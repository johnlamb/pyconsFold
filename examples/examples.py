#!/usr/bin/env python3
"""Tests to make sure pyconsFold is installed properly"""
import os
import sys
import pyconsFold
from collections import namedtuple, defaultdict
if not 'Bio' in sys.modules:
    print("Utils support not available without extra packages (BioPython, scipy and numpy")
    print("Install them manually or by using `pip3 install -U pyconsFold[all]`")
    sys.exit()
from pyconsFold.utils import npz_to_casp, pdb_to_npz, pdb_to_casp

### Setup for bin_values
Bin_values = namedtuple("Bin_values", ["bin_step", "min_bin_value", "max_bin_value", "ending"])
my_bin_values = {"dist":  Bin_values(0.5, 2, 16, ".rr"),         # In Angstrom
                  "omega": Bin_values(15, -180, 180, ".omega"),  # Dihedral angle in degrees
                  "theta": Bin_values(15, -180, 180, ".theta"),  # Dihedral angle in degrees
                  "phi":   Bin_values(15,    0, 180, ".phi")}    # Planar angle in degrees

###########################################################
### Convert npz to casp(rr format), use your own bin values
###########################################################
npz_to_casp("input_npz", fasta_file="fasta_file")  # Using default values
npz_to_casp("input_npz", fasta_file="fasta_file", info=["dist", "omega"], bin_values=my_bin_values, min_sep=5)  # Only convert distances and omega angles, 
                                                                                                      # Use custom bin_values and separation

#####################################
### Convert pdb/mmCIF into npz-format
#####################################
def pdb_to_npz(npz_name, pdb_file=False, mmCIF_file=False, std=1,
               DIST_STEP=0.5, OMEGA_STEP=15, THETA_STEP=15, PHI_STEP=15):
pdb_to_npz("outfile_name.npz", pdb_file="structure.pdb")  # Using defaults
pdb_to_npz("outfile_name.npz", pdb_file="structure.pdb", DIST_STEP=1, THETA_STEP=30, std=0.8)  # Using specific step for distance and theta and standard deviation per bin

#####################################
### Convert pdb/mmCIF into casp-format
#####################################
pdb_to_casp("outfile_name.rr", pdb_file="structure.pdb")  # Using defaults
pdb_to_casp("outfile_name.rr", pdb_file="structure.pdb", cutoff=8, confidence=0.8, sd=0.8)  # Using specific cutoff, confidence and standard deviation per bin
