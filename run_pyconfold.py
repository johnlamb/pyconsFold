#!/usr/bin/env python3
"""Docstring"""
import sys
import pyconfold

fa_file = sys.argv[1]
ss_file = sys.argv[2]
rr_file = sys.argv[3]
out_folder = sys.argv[4]
dist = True if sys.argv[5] == "True" else False

pyconfold.pyconfold(fa_file, ss_file, rr_file, out_folder, dist=dist, selectrr="1.5L")
