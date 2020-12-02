#!/usr/bin/env python3
"""Tests to make sure pyconsFold is installed properly"""
import os
import pyconsFold

passed = 1
total = 5
pyconsFold.model(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/short.fasta"),
                 os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/short.rr"),
                 "out_tests/modelout",
                 ss=os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/short.ss"), rr_sep=5)
print("Passed test {} of {}".format(passed, total))
passed += 1
pyconsFold.model_dist(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.rr"),
                      "out_tests/model_distout",)
print("Passed test {} of {}".format(passed, total))
passed += 1

pyconsFold.model_dist(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.npz"),
                      "out_tests/model_dist_npzout")
print("Passed test {} of {}".format(passed, total))
passed += 1

pyconsFold.model_dist(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.npz"),
                      "out_tests/model_dist_npz_anglesout", use_angles=True)
print("Passed test {} of {}".format(passed, total))
passed += 1

pyconsFold.model_dock(os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUA.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIUB.fasta"),
                      os.path.join(os.path.dirname(os.path.realpath(__file__)),"data/2HIU_dock.rr"),
                      "out_tests/model_dock")
print("Passed test {} of {}".format(passed, total))
