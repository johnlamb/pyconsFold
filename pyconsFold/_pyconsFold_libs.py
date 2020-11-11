#!/usr/bin/env python3
import os
from ._pyconsFold_helpers import print2file


def get_dist_neg_pos(mean, devi, lbd, lbd_dev=None):
    if lbd_dev is None:
        lbd_dev = round(lbd * devi, 1)
    return "{} {} {}".format(str(mean), str(lbd_dev), str(lbd_dev))


def load_ss_restraints(lbd, log_reference):
    # T      Helix or Parallel or anti-parallel or Unknown Strand Type
    # A1_A2  Atom1-Atom2 pair
    # Ref    Hydrogen bonding connector atom (reference hbond connector)
    # N      Neighborhood residue shifting on the hbond connector of R2 side.
    # For example, If R1:N and R2:O have hbond and S = +1, the restraint A1-A2
    # are for R1 and (R2+1)
    # Note: hbond distances are the distances between Nitrogen and Oxygen
    # Places to verify:
    # http://www.beta-sheet.org/page29/page51/page53/index.html
    # In this model, the HP sheet is composed of identical straight helical
    # chains with phi = -122 degrees, psi = 135 degrees, and a slightly
    # non-linear interchain H-bond angle delta of 170 degrees.
    # http://en.wikipedia.org/wiki/Alpha_helix
    # Residues in α-helices typically adopt backbone (φ, ψ) dihedral angles
    # around (-60°, -45°), as shown in the image at right
    # the H to O distance is about 2 Å (0.20 nm

    # T A Mean Standard_deviation
    res_dihe = {}
    res_dihe["A PSI"] = "136.91 " + str(round(lbd * 17.39, 3))
    res_dihe["A PHI"] = "-120.89 " + str(round(lbd * 21.98, 3))
    res_dihe["P PSI"] = "130.96 " + str(round(lbd * 16.66, 3))
    res_dihe["P PHI"] = "-115.00 " + str(round(lbd * 20.31, 3))
    res_dihe["U PSI"] = "134.95 " + str(round(lbd * 17.65, 3))
    res_dihe["U PHI"] = "-118.91 " + str(round(lbd * 21.73, 3))
    res_dihe["H PSI"] = "-41.51 " + str(round(lbd * 9.84, 3))
    res_dihe["H PHI"] = "-63.47 " + str(round(lbd * 9.20, 3))

    # T A1_A2 Ref N Mean Standard_deviation
    res_dist = {}
    res_dist["A O-O O +1"] = get_dist_neg_pos(7.73, 0.59, lbd)
    res_dist["A O-O O -1"] = get_dist_neg_pos(4.84, 0.16, lbd)
    res_dist["A O-O O 0"] = get_dist_neg_pos(3.57, 0.28, lbd)
    res_dist["A O-O H +1"] = get_dist_neg_pos(7.76, 0.60, lbd)
    res_dist["A O-O H -1"] = get_dist_neg_pos(4.90, 0.45, lbd)
    res_dist["A O-O H 0"] = get_dist_neg_pos(3.58, 0.31, lbd)
    res_dist["A C-C O +1"] = get_dist_neg_pos(7.66, 0.52, lbd)
    res_dist["A C-C O -1"] = get_dist_neg_pos(4.80, 0.17, lbd)
    res_dist["A C-C O 0"] = get_dist_neg_pos(4.96, 0.21, lbd)
    res_dist["A C-C H +1"] = get_dist_neg_pos(7.65, 0.51, lbd)
    res_dist["A C-C H -1"] = get_dist_neg_pos(4.85, 0.34, lbd)
    res_dist["A C-C H 0"] = get_dist_neg_pos(4.96, 0.21, lbd)
    res_dist["A N-N O +1"] = get_dist_neg_pos(5.09, 0.34, lbd)
    res_dist["A N-N O -1"] = get_dist_neg_pos(6.86, 0.40, lbd)
    res_dist["A N-N O 0"] = get_dist_neg_pos(4.42, 0.24, lbd)
    res_dist["A N-N H +1"] = get_dist_neg_pos(5.04, 0.21, lbd)
    res_dist["A N-N H -1"] = get_dist_neg_pos(6.85, 0.45, lbd)
    res_dist["A N-N H 0"] = get_dist_neg_pos(4.43, 0.25, lbd)
    res_dist["A CA-CA O +1"] = get_dist_neg_pos(6.43, 0.41, lbd)
    res_dist["A CA-CA O -1"] = get_dist_neg_pos(5.67, 0.28, lbd)
    res_dist["A CA-CA O 0"] = get_dist_neg_pos(5.26, 0.24, lbd)
    res_dist["A CA-CA H +1"] = get_dist_neg_pos(6.38, 0.36, lbd)
    res_dist["A CA-CA H -1"] = get_dist_neg_pos(5.71, 0.40, lbd)
    res_dist["A CA-CA H 0"] = get_dist_neg_pos(5.27, 0.25, lbd)
    res_dist["P O-O O +1"] = get_dist_neg_pos(7.90, 0.61, lbd)
    res_dist["P O-O O -1"] = get_dist_neg_pos(4.86, 0.16, lbd)
    res_dist["P O-O O 0"] = get_dist_neg_pos(3.78, 0.34, lbd)
    res_dist["P O-O H +1"] = get_dist_neg_pos(4.92, 0.40, lbd)
    res_dist["P O-O H -1"] = get_dist_neg_pos(8.02, 0.60, lbd)
    res_dist["P O-O H 0"] = get_dist_neg_pos(3.78, 0.32, lbd)
    res_dist["P C-C O +1"] = get_dist_neg_pos(8.03, 0.51, lbd)
    res_dist["P C-C O -1"] = get_dist_neg_pos(4.82, 0.17, lbd)
    res_dist["P C-C O 0"] = get_dist_neg_pos(5.21, 0.25, lbd)
    res_dist["P C-C H +1"] = get_dist_neg_pos(4.88, 0.34, lbd)
    res_dist["P C-C H -1"] = get_dist_neg_pos(7.87, 0.44, lbd)
    res_dist["P C-C H 0"] = get_dist_neg_pos(5.22, 0.22, lbd)
    res_dist["P N-N O +1"] = get_dist_neg_pos(8.14, 0.35, lbd)
    res_dist["P N-N O -1"] = get_dist_neg_pos(4.86, 0.40, lbd)
    res_dist["P N-N O 0"] = get_dist_neg_pos(5.13, 0.32, lbd)
    res_dist["P N-N H +1"] = get_dist_neg_pos(4.80, 0.18, lbd)
    res_dist["P N-N H -1"] = get_dist_neg_pos(7.54, 0.69, lbd)
    res_dist["P N-N H 0"] = get_dist_neg_pos(5.10, 0.28, lbd)
    res_dist["P CA-CA O +1"] = get_dist_neg_pos(8.55, 0.37, lbd)
    res_dist["P CA-CA O -1"] = get_dist_neg_pos(4.90, 0.29, lbd)
    res_dist["P CA-CA O 0"] = get_dist_neg_pos(6.21, 0.26, lbd)
    res_dist["P CA-CA H +1"] = get_dist_neg_pos(4.90, 0.28, lbd)
    res_dist["P CA-CA H -1"] = get_dist_neg_pos(7.49, 0.60, lbd)
    res_dist["P CA-CA H 0"] = get_dist_neg_pos(6.24, 0.24, lbd)
    res_dist["H O-O O +1"] = get_dist_neg_pos(8.40, 0.27, lbd)
    res_dist["H O-O O -1"] = get_dist_neg_pos(4.99, 0.16, lbd)
    res_dist["H O-O O 0"] = get_dist_neg_pos(6.12, 0.26, lbd)
    res_dist["H O-O H +1"] = get_dist_neg_pos(5.03, 0.31, lbd)
    res_dist["H O-O H -1"] = get_dist_neg_pos(8.43, 0.32, lbd)
    res_dist["H O-O H 0"] = get_dist_neg_pos(6.12, 0.26, lbd)
    res_dist["H C-C O +1"] = get_dist_neg_pos(8.16, 0.24, lbd)
    res_dist["H C-C O -1"] = get_dist_neg_pos(4.87, 0.13, lbd)
    res_dist["H C-C O 0"] = get_dist_neg_pos(6.09, 0.23, lbd)
    res_dist["H C-C H +1"] = get_dist_neg_pos(4.89, 0.23, lbd)
    res_dist["H C-C H -1"] = get_dist_neg_pos(8.17, 0.25, lbd)
    res_dist["H C-C H 0"] = get_dist_neg_pos(6.09, 0.23, lbd)
    res_dist["H N-N O +1"] = get_dist_neg_pos(8.07, 0.23, lbd)
    res_dist["H N-N O -1"] = get_dist_neg_pos(4.84, 0.19, lbd)
    res_dist["H N-N O 0"] = get_dist_neg_pos(6.10, 0.20, lbd)
    res_dist["H N-N H +1"] = get_dist_neg_pos(4.81, 0.13, lbd)
    res_dist["H N-N H -1"] = get_dist_neg_pos(8.08, 0.21, lbd)
    res_dist["H N-N H 0"] = get_dist_neg_pos(6.10, 0.20, lbd)
    res_dist["H CA-CA O +1"] = get_dist_neg_pos(8.63, 0.28, lbd)
    res_dist["H CA-CA O -1"] = get_dist_neg_pos(5.13, 0.20, lbd)
    res_dist["H CA-CA O 0"] = get_dist_neg_pos(6.16, 0.26, lbd)
    res_dist["H CA-CA H +1"] = get_dist_neg_pos(5.14, 0.21, lbd)
    res_dist["H CA-CA H -1"] = get_dist_neg_pos(8.64, 0.26, lbd)
    res_dist["H CA-CA H 0"] = get_dist_neg_pos(6.16, 0.26, lbd)

    # T Mean Standard_deviation
    res_strnd_OO = {}
    res_strnd_OO["A"] = get_dist_neg_pos(4.57, 0.30, lbd)
    res_strnd_OO["P"] = get_dist_neg_pos(4.57, 0.29, lbd)
    res_strnd_OO["U"] = get_dist_neg_pos(4.57, 0.30, lbd)

    # T Mean Standard_deviation
    res_hbnd = {}
    res_hbnd["A"] = get_dist_neg_pos(2.92, 0.16, lbd)
    res_hbnd["P"] = get_dist_neg_pos(2.93, 0.16, lbd)
    res_hbnd["H"] = get_dist_neg_pos(2.99, 0.17, lbd)

    if os.path.isfile(log_reference):
        os.remove(log_reference)

    log_ref = []
    for key, value in res_dihe.items():
        log_ref.append(key + '\t' + value)
    for key, value in res_strnd_OO.items():
        log_ref.append(key + '\t' + value)
    for key, value in res_dist.items():
        log_ref.append(key + '\t' + value)
    for key, value in res_hbnd.items():
        log_ref.append(key + '\t' + value)
    print2file(log_reference, '\n'.join(log_ref))

    return res_dihe, res_strnd_OO, res_dist, res_hbnd
