#!/usr/bin/env python3
"""Docstring"""
import sys
import pyconsFold
# import argparse

# parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument("fa_file", type=str, help="input sequence file in FASTA format")
# parser.add_argument("contact_file", type=str, help="Contact file in CASP or trRosetta npz format")
# parser.add_argument("out_dir", type=str, help="Output directory to write results")
# 
# parser.add_argument("-fa2_file", "--fa2", type=str, help="input sequence file in FASTA format for the second sequence")
# parser.add_argument("-ss_file", "--ss", default='', type=str, help="Secondary structure file in ss2/3 format")
# parser.add_argument("-L", "--L", default=0.55, type=float, help="Cutoff for contacts")
# 
# args = parser.parse_args()

########## Only uncomment ONE of these lines
pyconsFold.model("tests/short.fasta", "tests/short.rr", "out/modelout", ss="tests/short.ss", rr_sep=5, pcons=True)
pyconsFold.model_dist("tests/data/2HIUA.fasta", "tests/data/2HIUA.rr", "out/model_distout", pcons=True, tmscore_pdb_file="tests/data/2hiuA.pdb")
pyconsFold.model_dist("tests/data/2HIUA.fasta", "tests/data/2HIUA.npz", "out/model_dist_npzout", pcons=True, tmscore_pdb_file="tests/data/2hiuA.pdb")
pyconsFold.model_dist("tests/data/2HIUA.fasta", "tests/data/2HIUA.npz", "out/model_dist_npz_anglesout", use_angles=True, pcons=True, tmscore_pdb_file="tests/data/2hiuA.pdb")
pyconsFold.model_dock("tests/data/2HIUA.fasta", "tests/data/2HIUB.fasta", "tests/data/2HIU_dock.rr", "out/model_dock", pcons=True)

