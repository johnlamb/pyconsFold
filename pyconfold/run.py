#!/usr/bin/env python3
import argparse
import datetime
import glob
import os
import shutil
import sys

from ._pyconfold_helpers import (assess_dgsa, build_extended, build_models,
                                 check_programs, clean_output_dir,
                                 contact_restraints, process_arguments,
                                 sec_restraints, xyz_pdb, angle_restraints)
from ._pyconfold_libs import load_ss_restraints

from ._arguments import get_args

################## Default dssp ########################################
program_dssp = os.path.join(os.path.dirname(os.path.realpath(__file__)),\
                            "dssp/dssp-2.0.4-linux-amd64")
########## CNS from env variable #######################################
cns_suite = os.environ["CNS_SOLVE"]
cns_executable = cns_suite + "/intel-x86_64bit-linux/bin/cns_solve"
########################################################################
ATOMTYPE = {"CA": 1, "N": 1, "C": 1, "O": 1}
SHIFT = {"0": 1, "+1": 1, "-1": 1}


def initialize(fasta, rr, out_dir, dist, fasta2='', ss='', rrtype='cb', selectrr="all", omega='', theta='', mcount=20,
                      lbd=0.4, contwt=10, sswt=5, rep2=0.85, rr_pthres=0.0, rr_sep=5, debug=False):
    rep1 = 1.0
    atomselect = 2
    mode = "trial"

    base_dir = os.getcwd()
    dir_out = os.path.join(base_dir, out_dir)

    ### Check so that DSSP and CNS is accessable ###
    check_programs(program_dssp, cns_suite, cns_executable)

    ### Process all arguments, check validity and set up inputs ###
    fasta_file, fasta2_file, rr_file, ss_file, omega_file, theta_file, residues, f_id, selectrr, mini =\
            process_arguments(fasta, rr, out_dir, ss, rrtype, selectrr, fasta2, omega, theta, mcount, lbd, contwt, sswt, rep2, rr_pthres, rr_sep, debug)
    return fasta_file, fasta2_file, rr_file, dir_out, ss_file, omega_file, theta_file, residues, lbd,\
            selectrr, rrtype, contwt, sswt, mcount, mode, rep1, rep2, mini, f_id, atomselect


def setup_working_tree(dir_out):
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


def gen_initial_model(fasta_file, residues, fasta2_file, debug):
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
def model(fasta, rr, out_dir, fasta2='', ss='', dist=False, rr_pthres=0.0, top_models=5, save_step=False, debug=False):
    ############# 1. Setup variables and input files #########################
    fasta_file, fasta2_file, rr_file, out_dir, ss_file, omega_file,\
            theta_file, residues, lbd, selectrr, rrtype,\
            contwt, sswt, mcount, mode, rep1, rep2, mini, f_id, atomselect = initialize(fasta, rr, out_dir, dist, fasta2=fasta2, rr_pthres=rr_pthres, debug=debug)
    setup_working_tree(out_dir)

    ### Change working directory into stage1 for further work ###
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
    gen_initial_model(fasta_file, residues, fasta2_file, debug)

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
        ### Use angles files ###
        use_omega = True if omega_file else False
        use_theta = True if theta_file else False
        if use_omega or use_theta:
            angle_restraints(omega_file, use_omega, theta_file, use_theta, residues)

        ### Use the generated restraint files to build models ###
        build_models(stage, fasta_file, ss_file, contwt, sswt, mcount,
                     mode, rep1, rep2, mini, f_id, atomselect, out_dir, cns_suite, debug)

        ### Assess the models and rank them based on NOE energy ###
        assess_dgsa(stage, fasta_file, ss_file, out_dir, mcount, f_id, top_models,
                    program_dssp, debug)

        ##################### 3. Cleanup ##########################
    if not save_step:
        clean_output_dir(out_dir)

    if debug:
        print("\nFinished [{}]: {}".format("Pyconfold",
                                           datetime.datetime.now().
                                           replace(microsecond=0)))


### Model using distances, distance flag and rr_pthres is set ###
def model_dist(fasta, rr, out_dir, ss='', rr_pthres=0.55, top_models=5, save_step=False, debug=False):
    model(fasta, rr, out_dir, ss, dist=True, rr_pthres=rr_pthres, top_models=top_models, save_step=save_step, debug=debug)


### Docking, a second fasta file is supplied and distance is used by default (trRosetta output)
def model_dock(fasta, fasta2, rr, out_dir, dist=True, ss='', top_models=5, save_step=False, debug=False):
    model(fasta, rr, out_dir, ss, fasta2=fasta2, top_models=top_models, save_step=save_step, debug=debug)
    
