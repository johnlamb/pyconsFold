#!/usr/bin/env python3
import argparse
import datetime
import glob
import os
import shutil
import sys

from ._pyconsFold_helpers import (assess_dgsa, build_extended, build_models,
                                 check_programs, clean_output_dir,
                                 contact_restraints, process_arguments,
                                 sec_restraints, xyz_pdb, angle_restraints,
                                 qa_pcons, qa_tmscore, print2file)
from ._pyconsFold_libs import load_ss_restraints
from ._arguments import get_args

################## Default dssp ########################################
raw_program_dssp = os.path.join(os.path.dirname(os.path.realpath(__file__)),\
                            "dssp/dssp-2.0.4-linux-amd64")
raw_program_pcons = os.path.join(os.path.dirname(os.path.realpath(__file__)),\
                            "pcons/pcons")
raw_program_tmscore = os.path.join(os.path.dirname(os.path.realpath(__file__)),\
                            "TMscore/TMscore")
########## CNS from env variable #######################################
raw_cns_suite = os.environ["CNS_SOLVE"]
raw_cns_executable = raw_cns_suite + "/intel-x86_64bit-linux/bin/cns_solve"
########################################################################
### Check so that DSSP and CNS is accessable ###
program_dssp, program_pcons, program_tmscore, cns_suite, cns_executable = check_programs(raw_program_dssp, raw_program_pcons, raw_program_tmscore, raw_cns_suite, raw_cns_executable)

ATOMTYPE = {"CA": 1, "N": 1, "C": 1, "O": 1}
SHIFT = {"0": 1, "+1": 1, "-1": 1}


def _initialize(fasta, out_dir, dist, rr='', npz='', fasta2='', ss='', rrtype='cb', selectrr="all", omega='', theta='', mcount=20,
                      lbd=0.4, contwt=10, sswt=5, rep2=0.85, rr_pthres=0.0, rr_sep=5, rep1=1, atomselect=2, mode='trial',debug=False, bin_values={},tmscore_pdb_file=False):

    #### If we should use tmscore to compare to native structure, make sure it is installed and accessable ###
    base_dir = os.getcwd()
    dir_out = os.path.join(base_dir, out_dir)

    ### Process all arguments, check validity and set up inputs ###
    fasta_file, fasta2_file, rr_file, ss_file, omega_file, theta_file, residues, f_id, selectrr, mini =\
            process_arguments(fasta, rr, npz, dir_out, ss, rrtype, selectrr, fasta2, omega, theta, mcount, lbd, contwt, sswt, rep2, rr_pthres, rr_sep, debug)
    return fasta_file, fasta2_file, rr_file, dir_out, ss_file, omega_file, theta_file, residues, lbd,\
            selectrr, rrtype, contwt, sswt, mcount, mode, rep1, rep2, mini, f_id, atomselect, debug, tmscore_pdb_file


def _setup_working_tree(dir_out):
    ### Setup the work tree, copy from input ###
    module_base = os.path.dirname(os.path.realpath(__file__))

    shutil.copytree(os.path.join(dir_out, "input"),
                    os.path.join(dir_out, "stage1"))
    shutil.copy(os.path.join(module_base, "templates/gseq.inp"),
                os.path.join(dir_out, "stage1/"))
    shutil.copy(os.path.join(module_base, "templates/extn.inp"),
                os.path.join(dir_out, "stage1/"))
    shutil.copy(os.path.join(module_base, "templates/dgsa.inp"),
                os.path.join(dir_out, "stage1/"))
    shutil.copy(os.path.join(module_base, "templates/scalecoolsetupedited"),
                os.path.join(dir_out, "stage1/"))
    shutil.copy(os.path.join(module_base, "templates/scalehotedited"),
                os.path.join(dir_out, "stage1/"))


def _gen_initial_model(fasta_file, residues, fasta2_file, debug):
    ### Generate initial mtf and pdb ###
    if debug:
        print("\nBuild extended mtf and pdb...")

    if not os.path.isfile("extended.pdb"):
            build_extended(fasta_file, cns_suite, cns_executable, fasta2_file)
    defined_atoms = xyz_pdb("extended.pdb", "all")

    ### Check everything was generated ok ###
    # Hydrogen or Nitrogen must be present for all atoms, in order to apply
    # hydrogen bond restraints
    for ix, res in residues.items():
        if str(ix) + " H" in defined_atoms or\
           str(ix) + " N" in defined_atoms:
            continue
        else:
            print("ERROR!! Something went wrong with residue {} ".format(ix) +
                  "extended must have H or N for each residue!")
            sys.exit()

### Main model function, all other is based on this ###
def _model(fasta, out_dir, rr='', npz='', fasta2='', dist=False, rr_pthres=0.55, top_models=5, save_step=False, **kwargs):
    ############# 1. Setup variables and input files #########################
    fasta_file, fasta2_file, rr_file, out_dir, ss_file, omega_file,\
            theta_file, residues, lbd, selectrr, rrtype,\
            contwt, sswt, mcount, mode, rep1, rep2, mini, f_id, atomselect, debug, tmscore_pdb_file = _initialize(fasta, out_dir, dist, rr=rr, npz=npz, fasta2=fasta2, rr_pthres=rr_pthres, **kwargs)
    _setup_working_tree(out_dir)

    ### Change working directory into stage1 for further work ###
    start_dir = os.getcwd()
    os.chdir(os.path.join(out_dir, "stage1"))

    ### Only get default ss geometry data if we are going to use it ###
    if ss_file:
        res_dihe, res_strnd_OO, res_dist, res_hbnd =\
                load_ss_restraints(lbd, "ssrestraints.log")

    ############# 2. Start of Modeling, generate inital files and model #######
    if debug:
        print("\nStart [{}]: {}".format("Pyconfold",
                                        datetime.datetime.now().
                                        replace(microsecond=0)))
    _gen_initial_model(fasta_file, residues, fasta2_file, debug)

    ### Currently, only stage1 is implemented ###
    stage_list = []
    stage_list = ["stage1"]

    for stage in stage_list:
        os.chdir(os.path.join(out_dir, stage))

        if debug:
            print("\nStart {} job...".format(stage))

        contact_restraints(stage, selectrr, rrtype, out_dir, dist, rr_file)

        ### Use secondary structure if supplied ###
        if ss_file:
            sec_restraints(stage, ss_file, res_dihe, res_hbnd, res_dist,
                           res_strnd_OO, residues, ATOMTYPE, SHIFT, debug)
            ### Update the contact.tbl file with new number of restraints
            extra_restraints = 0
            for tbl in ["hbond.tbl", "dihedral.tbl", "ssnoe.tbl"]:
                if os.path.isfile(tbl):
                    with open(tbl,'r') as tbl_handle:
                        num = len(tbl_handle.read().split('\n'))
                        extra_restraints += num
            contact_to_write = ''
            with open("contact.tbl", 'r') as contact_handle:
                curr_restraints = int(contact_handle.readline().split('=')[1])
                contact_to_write += 'nrestraints=' + str(extra_restraints + curr_restraints) + '\n'
                for line in contact_handle:
                    contact_to_write += line
            os.remove("contact.tbl")
            print2file("contact.tbl", contact_to_write)

        ### Use angles files ###
        use_omega = True if omega_file else False
        use_theta = True if theta_file else False
        if use_omega or use_theta:
            angle_restraints(omega_file, use_omega, theta_file, use_theta, residues)

        ### Use the generated restraint files to build models ###
        build_models(stage, fasta_file, fasta2_file, ss_file, contwt, sswt, mcount,
                     mode, rep1, rep2, mini, f_id, atomselect, out_dir, cns_suite, debug)

        ### Assess the models and rank them based on NOE energy ###
        noe_scores = assess_dgsa(stage, fasta_file, ss_file, out_dir, mcount, f_id, top_models,
                    program_dssp, debug)
        ### Go back to star dir and run QA ###
        os.chdir(start_dir)
        ### Run QA programs for the top models ###
        pcons_scores = qa_pcons(out_dir, f_id, program_dssp, debug)

        ### If we have native structure, run with TMscore ###
        if tmscore_pdb_file:
            tmscores = qa_tmscore(out_dir, tmscore_pdb_file, debug)

        scores = ""
        for noe in noe_scores:
            scores += "NOE {} {}\n".format(noe[0], noe[1])
        for pcons in pcons_scores:
            scores += "pcons {} {}\n".format(pcons[0], pcons[1])
        if tmscore_pdb_file:
            for tmscore in tmscores:
                scores += "TMscore {} {}\n".format(tmscore[0], tmscore[1])

        print2file(out_dir + "/qa_scores.txt", scores)
        ##################### 3. Cleanup ##########################
    if not save_step:
        clean_output_dir(out_dir)

    if debug:
        print("\nFinished [{}]: {}".format("Pyconfold",
                                           datetime.datetime.now().
                                           replace(microsecond=0)))


### Base model function, uses simple binary contacts ###
def model(fasta, out_dir, rr='', npz='', top_models=5, save_step=False, **kwargs):
    _model(fasta, out_dir, rr=rr, npz=npz, dist=False,top_models=top_models, save_step=save_step, **kwargs)

### Model using distances, distance flag and rr_pthres is set ###
def model_dist(fasta, out_dir, rr='', npz='', rr_pthres=0.55, top_models=5, save_step=False, debug=False, **kwargs):
    _model(fasta, out_dir, rr=rr, npz=npz, dist=True, rr_pthres=rr_pthres, top_models=top_models, save_step=save_step, **kwargs)


### Docking, a second fasta file is supplied and distance is used by default (trRosetta output)
def model_dock(fasta, fasta2, out_dir, rr='', npz='', dist=True, top_models=5, save_step=False, **kwargs):
    _model(fasta, out_dir, rr=rr, npz=npz, fasta2=fasta2, top_models=top_models, save_step=save_step, **kwargs)
    
