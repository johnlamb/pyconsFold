#!/usr/bin/env python3
"""Docstring"""
import sys
import pyconfold

fa_file = sys.argv[1]
rr_file = sys.argv[2]
ss_file = sys.argv[3]
out_folder = sys.argv[4]
# dist = False
# dist_error = False
# if len(sys.argv) > 5:
#     dist = True if sys.argv[5] == "True" else False
#     dist_error = True if sys.argv[6] == "True" else False

# pyconfold.pyconfold(debug=False)
pyconfold.model(fa_file, rr_file, out_folder, ss=ss_file)
