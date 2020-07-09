#!/usr/bin/env python3
"""Docstring"""
import sys
import pyconfold
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("fa_file", type=str, help="input sequence file in FASTA format")
parser.add_argument("rr_file", type=str, help="Contact file in CASP format, with or without error as last column")
parser.add_argument("out_dir", type=str, help="Output directory to write results")

parser.add_argument("-fa2_file", "--fa2", type=str, help="input sequence file in FASTA format for the second sequence")
parser.add_argument("-ss_file", "--ss", type=str, help="Secondary structure file in ss2/3 format")

args = parser.parse_args()
# dist = False
# dist_error = False
# if len(sys.argv) > 5:
#     dist = True if sys.argv[5] == "True" else False
#     dist_error = True if sys.argv[6] == "True" else False

########## Only uncomment ONE of these lines
# pyconfold.model(args.fa_file, args.rr_file, args.out_folder, ss=args.ss_file)
pyconfold.model_dist(args.fa_file, args.rr_file, args.out_folder)
# pyconfold.model_dock(args.fa_file, args.fa2_file, args.rr_file, args.out_folder)
