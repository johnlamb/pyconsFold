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
program_dssp = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dssp/dssp-2.0.4-linux-amd64")
########## The following three for Keb #################################
cns_suite = os.environ["CNS_SOLVE"]
cns_executable = cns_suite + "/intel-x86_64bit-linux/bin/cns_solve"
########################################################################
# pair=None  not implemented yet

# def pyconfold(fasta, ss, rr, dir_out, save_steps=False, num_top_models=5, rrtype="cb", mcount=20, selectrr="all",
#               lbd=0.4, contwt=10, sswt=5, rep2=0.85, pthres=7.0, dist=False, dist_error=False, debug=False):
def pyconfold(debug=False):
    args = get_args()

    rep1 = 1.0
    mini = 0
    # pthres_hbond = args.pthres + 3.0
    pthres_hbond = args.pthres + 3.0
    atomselect = 2
    mode = "trial"

    ATOMTYPE = {"CA": 1, "N": 1, "C": 1, "O": 1}
    SHIFT = {"0": 1, "+1": 1, "-1": 1}

    base_dir = os.getcwd()
    # dir_out = os.path.join(base_dir, args.dir_out)
    dir_out = os.path.join(base_dir, args.out_dir)

    check_programs(program_dssp, cns_suite, cns_executable)
    # fasta_file, rr_file, ss_file, pair_file, residues, f_id, selectrr, mini =\
    #         process_parameters(args)
    # pair  not implemented
    if not os.path.isfile(args.fasta):
        print("ERROR! Fasta file {} does not exist!".format(args.fasta))
        sys.exit()
    if not os.path.isfile(args.sec_struct):
        print("ERROR! Secondary structure file {} " +
              "does not exist!".format(args.sec_struct))
        sys.exit()
    if not os.path.isfile(args.rr):
        print("ERROR! Contact file {} does not exist!".format(args.rr))
        sys.exit()

    fasta_file, rr_file, ss_file, omega_file, theta_file, residues, f_id, selectrr, mini =\
            process_arguments(args.fasta, args.sec_struct, args.rr, args.out_dir, args.rr_type, args.omega, args.theta, args.model_count, args.select_rr,
                              args.lbd, args.contwt, args.sswt, args.rep2, args.pthres, debug, args.rr_pthres, 5)
    # base_dir = os.path.dirname(os.path.realpath(__file__))
    # base_dir = os.getcwd()
    # print(base_dir)
    # dir_out = os.path.join(base_dir, args.dir_out)
    # module_base = os.path.join(base_dir, os.path.dirname(os.path.relpath(__file__)))
    module_base = os.path.dirname(os.path.realpath(__file__))
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    # print(base_dir)
    # print(dir_path)
    # print(module_base)
    # sys.exit()
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


    os.chdir(os.path.join(dir_out, "stage1"))

    # res_dihe, res_strnd_OO, res_dist, res_hbnd =\
    #         load_ss_restraints(args.lbd, "ssrestraints.log")
    res_dihe, res_strnd_OO, res_dist, res_hbnd =\
            load_ss_restraints(args.lbd, "ssrestraints.log")

    # print("\nStart [{}]: {}".format(parser.prog,
    #                                 datetime.datetime.now().
    #                                 replace(microsecond=0)))
    if debug:
        print("\nStart [{}]: {}".format("Pyconfold",
                                        datetime.datetime.now().
                                        replace(microsecond=0)))
        print("\nBuild extended mtf and pdb...")

    if not os.path.isfile("extended.pdb"):
            build_extended(fasta_file, cns_suite, cns_executable)
    defined_atoms = xyz_pdb("extended.pdb", "all")

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

    stage_list = []
    # if args.stage2:
    #     stage_list = ["stage1", "stage2"]
    # else:
    # Stage 2 is not implemented
    stage_list = ["stage1"]

    for stage in stage_list:
        # Stage 2 is not implemented
        # if stage == "stage2":
        #     os.mkdir(os.path.join(dir_out, stage))

        os.chdir(os.path.join(dir_out, stage))

        if debug:
            print("\nStart {} job...".format(stage))
        # Stage 2 is not implemented
        # if stage == "stage2":
        #     sys.exit()
        #     for filename in glob.glob(os.path.join(dir_out, "stage1",
        #                                                     "extended.*")):
        #         shutil.copy(filename, ".")
        #     shutil.copy(os.path.join(dir_out, "stage1", "{}.fasta".format(f_id)),
        #                 ".")
        #     shutil.copy(os.path.join(dir_out, "stage1", "{}.ss".format(f_id)),
        #                 ".")
        # contact_restraints(stage, selectrr, args.rrtype, rr_file)
        contact_restraints(stage, selectrr, args.rr_type, dir_out, False, rr_file)
        # if not args.noss:
        if ss_file:
            sec_restraints(stage, ss_file, res_dihe, res_hbnd, res_dist,
                           res_strnd_OO, residues, ATOMTYPE, SHIFT, debug)
        ###############################################################
        # angle_restraints(omega_file, args.use_omega, theta_file, args.use_theta, residues)
        ###############################################################
        # build_models(stage, fasta_file, ss_file, args.contwt, args.sswt,
        #              args.mcount, mode, rep1, args.rep2, mini, f_id, atomselect,
        #              dir_out, cns_suite)
        build_models(stage, fasta_file, ss_file, args.contwt, args.sswt,
                     args.model_count, mode, rep1, args.rep2, mini, f_id, atomselect,
                     dir_out, cns_suite, debug)
        # assess_dgsa(stage, fasta_file, ss_file, dir_out, args.mcount, f_id,
        #             program_dssp)
        assess_dgsa(stage, fasta_file, ss_file, dir_out, args.model_count, f_id, args.top_models,
                    program_dssp, debug)

    if not args.sstep:
        clean_output_dir(dir_out)

    # print("\nFinished [{}]: {}".format(parser.prog,
    #                                    datetime.datetime.now().
    #                                    replace(microsecond=0)))
    print("\nFinished [{}]: {}".format("Pyconfold",
                                       datetime.datetime.now().
                                       replace(microsecond=0)))


def pyconfold_dist(debug=False):
    args = get_args()

    rep1 = 1.0
    mini = 0
    pthres_hbond = args.pthres + 3.0
    atomselect = 2
    mode = "trial"

    ATOMTYPE = {"CA": 1, "N": 1, "C": 1, "O": 1}
    SHIFT = {"0": 1, "+1": 1, "-1": 1}

    base_dir = os.getcwd()
    # dir_out = os.path.join(base_dir, args.dir_out)
    dir_out = os.path.join(base_dir, args.out_dir)

    check_programs(program_dssp, cns_suite, cns_executable)
    # fasta_file, rr_file, ss_file, pair_file, residues, f_id, selectrr, mini =\
    #         process_parameters(args)
    # pair  not implemented
    if not os.path.isfile(args.fasta):
        print("ERROR! Fasta file {} does not exist!".format(args.fasta))
        sys.exit()
    if args.sec_struct and not os.path.isfile(args.sec_struct):
        print("ERROR! Secondary structure file {} " +
              "does not exist!".format(args.sec_struct))
        sys.exit()
    if not os.path.isfile(args.rr):
        print("ERROR! Contact file {} does not exist!".format(args.rr))
        sys.exit()

    fasta_file, rr_file, ss_file, omega_file, theta_file, residues, f_id, selectrr, mini =\
            process_arguments(args.fasta, args.sec_struct, args.rr, args.out_dir, args.rr_type, args.omega, args.theta, args.model_count, args.select_rr,
                              args.lbd, args.contwt, args.sswt, args.rep2, args.pthres, debug, args.rr_pthres, 0)
    # base_dir = os.path.dirname(os.path.realpath(__file__))
    # base_dir = os.getcwd()
    # print(base_dir)
    # dir_out = os.path.join(base_dir, args.dir_out)
    # module_base = os.path.join(base_dir, os.path.dirname(os.path.relpath(__file__)))
    module_base = os.path.dirname(os.path.realpath(__file__))
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    # print(base_dir)
    # print(dir_path)
    # print(module_base)
    # sys.exit()
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


    os.chdir(os.path.join(dir_out, "stage1"))

    # res_dihe, res_strnd_OO, res_dist, res_hbnd =\
    #         load_ss_restraints(args.lbd, "ssrestraints.log")
    res_dihe, res_strnd_OO, res_dist, res_hbnd =\
            load_ss_restraints(args.lbd, "ssrestraints.log")

    # print("\nStart [{}]: {}".format(parser.prog,
    #                                 datetime.datetime.now().
    #                                 replace(microsecond=0)))
    if debug:
        print("\nStart [{}]: {}".format("Pyconfold",
                                        datetime.datetime.now().
                                        replace(microsecond=0)))
        print("\nBuild extended mtf and pdb...")

    if not os.path.isfile("extended.pdb"):
            build_extended(fasta_file, cns_suite, cns_executable)
    defined_atoms = xyz_pdb("extended.pdb", "all")

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

    stage_list = []
    # if args.stage2:
    #     stage_list = ["stage1", "stage2"]
    # else:
    # Stage 2 is not implemented
    stage_list = ["stage1"]

    for stage in stage_list:
        # Stage 2 is not implemented
        # if stage == "stage2":
        #     os.mkdir(os.path.join(dir_out, stage))

        os.chdir(os.path.join(dir_out, stage))

        if debug:
            print("\nStart {} job...".format(stage))
        # Stage 2 is not implemented
        # if stage == "stage2":
        #     sys.exit()
        #     for filename in glob.glob(os.path.join(dir_out, "stage1",
        #                                                     "extended.*")):
        #         shutil.copy(filename, ".")
        #     shutil.copy(os.path.join(dir_out, "stage1", "{}.fasta".format(f_id)),
        #                 ".")
        #     shutil.copy(os.path.join(dir_out, "stage1", "{}.ss".format(f_id)),
        #                 ".")
        # contact_restraints(stage, selectrr, args.rrtype, rr_file)
        contact_restraints(stage, selectrr, args.rr_type, dir_out, True, rr_file)
        # if not args.noss:
        #     sec_restraints(stage, ss_file, res_dihe, res_hbnd, res_dist,
        #                    res_strnd_OO, residues, ATOMTYPE, SHIFT, debug)
        ###############################################################
        # angle_restraints(omega_file, args.use_omega, theta_file, args.use_theta, residues)
        ###############################################################
        # build_models(stage, fasta_file, ss_file, args.contwt, args.sswt,
        #              args.mcount, mode, rep1, args.rep2, mini, f_id, atomselect,
        #              dir_out, cns_suite)
        build_models(stage, fasta_file, ss_file, args.contwt, args.sswt,
                     args.model_count, mode, rep1, args.rep2, mini, f_id, atomselect,
                     dir_out, cns_suite, debug)
        # assess_dgsa(stage, fasta_file, ss_file, dir_out, args.mcount, f_id,
        #             program_dssp)
        assess_dgsa(stage, fasta_file, ss_file, dir_out, args.model_count, f_id, args.top_models,
                    program_dssp, debug)

    if not args.sstep:
        clean_output_dir(dir_out)

    # print("\nFinished [{}]: {}".format(parser.prog,
    #                                    datetime.datetime.now().
    #                                    replace(microsecond=0)))
    print("\nFinished [{}]: {}".format("Pyconfold",
                                       datetime.datetime.now().
                                       replace(microsecond=0)))


def pyconfold_dock(debug=False):
    args = get_args()

    rep1 = 1.0
    mini = 0
    pthres_hbond = args.pthres + 3.0
    atomselect = 2
    mode = "trial"

    ATOMTYPE = {"CA": 1, "N": 1, "C": 1, "O": 1}
    SHIFT = {"0": 1, "+1": 1, "-1": 1}

    base_dir = os.getcwd()
    # dir_out = os.path.join(base_dir, args.dir_out)
    dir_out = os.path.join(base_dir, args.out_dir)

    check_programs(program_dssp, cns_suite, cns_executable)
    # fasta_file, rr_file, ss_file, pair_file, residues, f_id, selectrr, mini =\
    #         process_parameters(args)
    # pair  not implemented
    if not os.path.isfile(args.fasta):
        print("ERROR! Fasta file {} does not exist!".format(args.fasta))
        sys.exit()
    if not os.path.isfile(args.fa2):
        print("ERROR! Fasta file #2 {} does not exist!".format(args.fa2))
        sys.exit()
    if args.sec_struct and not os.path.isfile(args.sec_struct):
        print("ERROR! Secondary structure file {} " +
              "does not exist!".format(args.sec_struct))
        sys.exit()
    if not os.path.isfile(args.rr):
        print("ERROR! Contact file {} does not exist!".format(args.rr))
        sys.exit()

    fasta_file, fasta2_file, rr_file, ss_file, omega_file, theta_file, residues, f_id, selectrr, mini =\
            process_arguments(args.fasta, args.sec_struct, args.rr, args.out_dir, args.rr_type, args.omega, args.theta, args.model_count, args.select_rr,
                              args.lbd, args.contwt, args.sswt, args.rep2, args.pthres, debug, args.rr_pthres, 0, args.fa2)
    # base_dir = os.path.dirname(os.path.realpath(__file__))
    # base_dir = os.getcwd()
    # print(base_dir)
    # dir_out = os.path.join(base_dir, args.dir_out)
    # module_base = os.path.join(base_dir, os.path.dirname(os.path.relpath(__file__)))
    module_base = os.path.dirname(os.path.realpath(__file__))
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    # print(base_dir)
    # print(dir_path)
    # print(module_base)
    # sys.exit()
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


    os.chdir(os.path.join(dir_out, "stage1"))

    # res_dihe, res_strnd_OO, res_dist, res_hbnd =\
    #         load_ss_restraints(args.lbd, "ssrestraints.log")
    res_dihe, res_strnd_OO, res_dist, res_hbnd =\
            load_ss_restraints(args.lbd, "ssrestraints.log")

    # print("\nStart [{}]: {}".format(parser.prog,
    #                                 datetime.datetime.now().
    #                                 replace(microsecond=0)))
    if debug:
        print("\nStart [{}]: {}".format("Pyconfold",
                                        datetime.datetime.now().
                                        replace(microsecond=0)))
        print("\nBuild extended mtf and pdb...")

    if not os.path.isfile("extended.pdb"):
            build_extended(fasta_file, cns_suite, cns_executable, fasta2_file)
    sys.exit()
    defined_atoms = xyz_pdb("extended.pdb", "all")

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

    stage_list = []
    # if args.stage2:
    #     stage_list = ["stage1", "stage2"]
    # else:
    # Stage 2 is not implemented
    stage_list = ["stage1"]

    for stage in stage_list:
        # Stage 2 is not implemented
        # if stage == "stage2":
        #     os.mkdir(os.path.join(dir_out, stage))

        os.chdir(os.path.join(dir_out, stage))

        if debug:
            print("\nStart {} job...".format(stage))
        # Stage 2 is not implemented
        # if stage == "stage2":
        #     sys.exit()
        #     for filename in glob.glob(os.path.join(dir_out, "stage1",
        #                                                     "extended.*")):
        #         shutil.copy(filename, ".")
        #     shutil.copy(os.path.join(dir_out, "stage1", "{}.fasta".format(f_id)),
        #                 ".")
        #     shutil.copy(os.path.join(dir_out, "stage1", "{}.ss".format(f_id)),
        #                 ".")
        # contact_restraints(stage, selectrr, args.rrtype, rr_file)
        contact_restraints(stage, selectrr, args.rr_type, dir_out, True, rr_file)
        # if not args.noss:
        #     sec_restraints(stage, ss_file, res_dihe, res_hbnd, res_dist,
        #                    res_strnd_OO, residues, ATOMTYPE, SHIFT, debug)
        ###############################################################
        # angle_restraints(omega_file, args.use_omega, theta_file, args.use_theta, residues)
        ###############################################################
        # build_models(stage, fasta_file, ss_file, args.contwt, args.sswt,
        #              args.mcount, mode, rep1, args.rep2, mini, f_id, atomselect,
        #              dir_out, cns_suite)
        build_models(stage, fasta_file, ss_file, args.contwt, args.sswt,
                     args.model_count, mode, rep1, args.rep2, mini, f_id, atomselect,
                     dir_out, cns_suite, debug)
        # assess_dgsa(stage, fasta_file, ss_file, dir_out, args.mcount, f_id,
        #             program_dssp)
        assess_dgsa(stage, fasta_file, ss_file, dir_out, args.model_count, f_id, args.top_models,
                    program_dssp, debug)

    if not args.sstep:
        clean_output_dir(dir_out)

    # print("\nFinished [{}]: {}".format(parser.prog,
    #                                    datetime.datetime.now().
    #                                    replace(microsecond=0)))
    print("\nFinished [{}]: {}".format("Pyconfold",
                                       datetime.datetime.now().
                                       replace(microsecond=0)))
