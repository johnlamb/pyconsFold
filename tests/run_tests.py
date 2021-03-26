#!/usr/bin/env python3
"""Tests to make sure pyconsFold is installed properly"""
import os
import sys
import pyconsFold

passed = 1
total = 6
pyconsFold.model(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/short.fasta"),
                 os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/short.rr"),
                 "out_tests/modelout_short",
                 ss=os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/short.ss"), rr_sep=5)
print("Passed test {} of {}".format(passed, total))
passed += 1
pyconsFold.model(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                 os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.rr"),
                 "out_tests/modelout",
                 cutoff=16, rr_sep=5)
print("Passed test {} of {}".format(passed, total))
passed += 1
pyconsFold.model_dist(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.rr"),
                      "out_tests/model_distout",)
print("Passed test {} of {}".format(passed, total))
passed += 1


if 'Bio' in sys.modules:
    pyconsFold.model_dist(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                          os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.npz"),
                          "out_tests/model_dist_npzout", dgsa_seed=9999)
    print("Passed test {} of {}".format(passed, total))
else:
    print("NPZ support not available without extra packages (BioPython, scipy and numpy")
    print("Install them manually or by using `pip3 install -U pyconsFold[all]`")
    print("Skipping test {}".format(passed + 1))
passed += 1

if 'Bio' in sys.modules:
    pyconsFold.model_dist(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                          os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.npz"),
                          "out_tests/model_dist_npz_anglesout", use_angles=True)
    print("Passed test {} of {}".format(passed, total))
else:
    print("NPZ support not available without extra packages (BioPython, scipy and numpy")
    print("Install them manually or by using `pip3 install -U pyconsFold[all]`")
    print("Skipping test {}".format(passed + 1))
passed += 1

pyconsFold.model_dock(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUB.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIU_dock.rr"),
                      "out_tests/model_dock")
print("Passed test {} of {}".format(passed, total))
