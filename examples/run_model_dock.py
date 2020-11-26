#!/usr/bin/env python3
"""Run docking modelling"""
import sys
import pyconsFold
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("fa_file", type=str, help="first input sequence file in FASTA format")
parser.add_argument("fa2_file", type=str, help="second input sequence file in FASTA format")
parser.add_argument("contact_file", type=str, help="Contact file in CASP or trRosetta npz format, combind contacts")
parser.add_argument("out_dir", type=str, help="Output directory to write results")
parser.add_argument("-stage2", "--s2", action="store_true", help="Also rung stage 2")

args = parser.parse_args()

pyconsFold.model_dock(args.fa_file, args.fa2_file, args.contact_file, args.out_dir, stage2=args.s2)

