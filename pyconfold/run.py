#!/usr/bin/env python3
import os
import argparse
import shutil
import glob
import sys
import datetime
from ._pyconfold_helpers import check_programs,\
                            process_parameters,\
                            process_arguments,\
                            build_extended, xyz_pdb,\
                            contact_restraints, sec_restraints, build_models,\
                            assess_dgsa, clean_output_dir
from ._pyconfold_libs import load_ss_restraints

# program_dssp = "/home/johnlamb/bin/dssp-2.0.4-linux-amd64"
program_dssp = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dssp/dssp-2.0.4-linux-amd64")

# pair=None  not implemented yet
def pyconfold(fasta, ss, rr, dir_out, save_steps=False, num_top_models=5, rrtype="cb", mcount=20, selectrr="all",
              lbd=0.4, contwt=10, sswt=5, rep2=0.85, pthres=7.0, debug=False):
    # cns_suite = "<enter CNS-path here>"
    cns_suite = os.environ["CNS_SOLVE"]
    cns_executable = cns_suite + "/intel-x86_64bit-linux/bin/cns_solve"

    # parser = argparse.ArgumentParser(description="Parse Confold arguments")
    # parser.add_argument("fasta", help="Input fasta file")
    # parser.add_argument("ss", help="Secondary structure file")
    # parser.add_argument("rr", help="Contact file")
    # parser.add_argument("dir_out", help="Output directory")
    # # parser.add_argument("-pair", help="Pair file")
    # parser.add_argument("-rrtype", default="cb", choices=["ca", "cb"],
    #                     help="Contact type")
    # # parser.add_argument("-stage2", default=0, type=int, choices=range(0, 4),
    # #                     help="Run stage 2")
    # parser.add_argument("-mcount", default=20, type=int, help="Model count")
    # parser.add_argument("-selectrr", default="all", help="Contacts to use")
    # parser.add_argument("-lbd", default=0.4, type=float, help="Lambda")
    # parser.add_argument("-contwt", default=10, type=float,
    #                     help="Contact restraint weights")
    # parser.add_argument("-sswt", default=5, type=float,
    #                     help="SS restraint weight")
    # parser.add_argument("-rep2", default=0.85, type=float, help="Rep2")
    # parser.add_argument("-pthres", default=7.0, type=float, help="Threshold")
    # args = parser.parse_args()

    rep1 = 1.0
    mini = 0
    # pthres_hbond = args.pthres + 3.0
    pthres_hbond = pthres + 3.0
    atomselect = 2
    mode = "trial"

    ATOMTYPE = {"CA": 1, "N": 1, "C": 1, "O": 1}
    SHIFT = {"0": 1, "+1": 1, "-1": 1}

    base_dir = os.getcwd()
    # dir_out = os.path.join(base_dir, args.dir_out)
    dir_out = os.path.join(base_dir, dir_out)

    check_programs(program_dssp, cns_suite, cns_executable)
    # fasta_file, rr_file, ss_file, pair_file, residues, f_id, selectrr, mini =\
    #         process_parameters(args)
    # pair  not implemented
    fasta_file, rr_file, ss_file, residues, f_id, selectrr, mini =\
            process_arguments(fasta, ss, rr, dir_out, rrtype, mcount, selectrr,
                              lbd, contwt, sswt, rep2, pthres, debug)
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
    shutil.copy(os.path.join(module_base,  "templates/gseq.inp"),
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
            load_ss_restraints(lbd, "ssrestraints.log")

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
    stage_list = ["stage1"]

    for stage in stage_list:
        if stage == "stage2":
            os.mkdir(os.path.join(dir_out, stage))

        os.chdir(os.path.join(dir_out, stage))

        print("\nStart {} job...".format(stage))
        if stage == "stage2":
            sys.exit()
            for filename in glob.glob(os.path.join(dir_out, "stage1",
                                                            "extended.*")):
                shutil.copy(filename, ".")
            shutil.copy(os.path.join(dir_out, "stage1", "{}.fasta".format(f_id)),
                        ".")
            shutil.copy(os.path.join(dir_out, "stage1", "{}.ss".format(f_id)),
                        ".")
        # contact_restraints(stage, selectrr, args.rrtype, rr_file)
        contact_restraints(stage, selectrr, rrtype, rr_file)
        sec_restraints(stage, ss_file, pair_file, res_dihe, res_hbnd, res_dist,
                       res_strnd_OO, residues, ATOMTYPE, SHIFT)
        # build_models(stage, fasta_file, ss_file, args.contwt, args.sswt,
        #              args.mcount, mode, rep1, args.rep2, mini, f_id, atomselect,
        #              dir_out, cns_suite)
        build_models(stage, fasta_file, ss_file, contwt, sswt,
                     mcount, mode, rep1, rep2, mini, f_id, atomselect,
                     dir_out, cns_suite)
        # assess_dgsa(stage, fasta_file, ss_file, dir_out, args.mcount, f_id,
        #             program_dssp)
        assess_dgsa(stage, fasta_file, ss_file, dir_out, mcount, f_id, num_top_models,
                    program_dssp)

    if not save_steps:
        clean_output_dir(dir_out)

    # print("\nFinished [{}]: {}".format(parser.prog,
    #                                    datetime.datetime.now().
    #                                    replace(microsecond=0)))
    print("\nFinished [{}]: {}".format("Pyconfold",
                                       datetime.datetime.now().
                                       replace(microsecond=0)))
