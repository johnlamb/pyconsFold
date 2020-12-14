import argparse
from Bio.PDB.vectors import rotaxis, calc_angle, calc_dihedral
from Bio.PDB.Polypeptide import is_aa
from math import pi
import numpy as np
import scipy.stats as st
from collections import namedtuple
import sys
import os


def npz_to_casp(input_file, info=["all"], fasta_file=None,  fasta2_file=None, out_base_path="",
        min_sep=0, pthres=0.15, bin_values={}):
    """Converts trRosetta distance and contacts file into CASP formated files

    Positional arguments:
    input_file    -- Input npz file

    Keyword arguments:
    info          -- What info to extract? A list of one or more of ['all']/['dist', 'omega', 'theta', 'phi']
    fasta_file    -- Fasta file to fill sequence information for the header 
    fasta2_file   -- Support for dual files in case of docking in pyconsFold
    out_base_path -- Output directory if other than current directory
    min_sep       -- Minimum separation between contacts
    pthres        -- Only extract contacts with a probability above this threshold
    bin_values    -- Named tuple of dictionary of bin values for the different infos, see below for defaults

    Bin value defaults:
    Bin_values = namedtuple("Bin_values", ["bin_step", "min_bin_value", "max_bin_value", "ending"])
    default_bin_values = {"dist":  Bin_values(0.5, 2, 16, ".rr"),    # in angstrom
                      "omega": Bin_values(15, -180, 180, ".omega"),  # dihedral angle in degrees
                      "theta": Bin_values(15, -180, 180, ".theta"),  # dihedral angle in degrees
                      "phi":   Bin_values(15,    0, 180, ".phi")}    # planar angle in degrees
    """

    Bin_values = namedtuple("Bin_values", ["bin_step", "min_bin_value", "max_bin_value", "ending"])
    default_bin_values = {"dist":  Bin_values(0.5, 2, 16, ".rr"),    # In Angstrom
                      "omega": Bin_values(15, -180, 180, ".omega"),  # Dihedral angle in degrees
                      "theta": Bin_values(15, -180, 180, ".theta"),  # Dihedral angle in degrees
                      "phi":   Bin_values(15,    0, 180, ".phi")}    # Planar angle in degrees

    ### Has the user supplied their own bin values?
    if bin_values:
        for key in bin_values.keys():
            if key in default_bin_values:
                default_bin_values[key] = Bin_values(*bin_values[key])
    if 'all' in info:
        info_keys = ['dist', 'omega', 'theta', 'phi']
    else:
        info_keys = info


    for info_key in info_keys:
        try:
            with np.load(input_file) as npz_file:
                raw_data = npz_file[info_key]
        except:
            print("Error reading npz file: " + input_file)
            print("Is file and info key {} correct?".format(info_key))
            sys.exit()
        target = '.'.join(input_file.split('/')[-1].split('.')[:-1])

        header = '\n'.join(("PFRMAT RR",
                           "TARGET {}".format(target),
                           "AUTHOR pyconsFold",
                           "METHOD trRosetta_contact",
                           "REMARK {}".format(info_key),
                           "MODEL 1\n"))
        # Length of the protein
        nres = raw_data.shape[0]

        # If a fasta file is supplied, add in the sequence to the header
        if fasta_file:
            n=50
            try:
                with open(fasta_file) as fasta_handle:
                    fasta_seq = ''.join(fasta_handle.read().split('\n')[1:]).strip()
            except:
                print("Error reading fasta file: "+fasta_file)
                sys.exit()

            if fasta2_file:
                try:
                    with open(fasta2_file) as fasta2_handle:
                        fasta_seq += ''.join(fasta2_handle.read().split('\n')[1:]).strip()
                except:
                    print("Error reading fasta file: "+fasta_file)
                    sys.exit()

            # Assert the fasta sequence is of the same length as the resulting contacts results
            fasta_length = len(fasta_seq)
            assert nres == fasta_length, "Sequence from fasta file is not of the same length as" +\
                                           " the predictions, {} vs {}".format(fasta_length, nres)

            header += '\n'.join(fasta_seq[i:i+n] for i in range(0, len(fasta_seq), n))
            header += '\n' 

        bin_values= default_bin_values[info_key] 

        # Number of bins, minus one for the first contact binary bin
        nbins = raw_data.shape[-1] - 1

        # # How many of the bins do we want?
        wanted_bins = int((bin_values.max_bin_value - bin_values.min_bin_value)/bin_values.bin_step)
        if wanted_bins > nbins:
            print("Number of wanted bins are higher than existing bins")
            sys.exit()
        # Generate up the midpoint of each distance bin to be able to
        # calculate correct distances. Measure in Angstrom
        bins = np.array([(bin_values.min_bin_value+bin_values.bin_step/2)+bin_values.bin_step*i for i in range(wanted_bins)])

        # Pull out the predictions, first array element is 
        # probability of non-contact, only pull out as many bins as we want
        # [L, L, nbins] -> [L, L, wanted_bins]
        raw_predictions = raw_data[:,:,1:wanted_bins+1]

        # What is the total probability of contact?
        contact_probabilities = np.sum(raw_predictions, axis=-1)

        # Normalize all probabilities by dividing by the sum of the raw_predictions
        # across the last dimension, [L, L, nbins] -> [L, L, nbins]
        norm_predictions = np.divide(raw_predictions, np.sum(raw_predictions, axis=-1)[:,:,None])

        # Multiply normalized predictions with the bin (sizes) to get the correct value distribution
        values = np.multiply(norm_predictions, bins[None, None, :])

        # Make some assertions to confirm array operations har occured correctly
        # Has the normalization worked?
        assert np.array_equal(np.divide(raw_predictions[0, 5, :],np.sum(raw_predictions[0,5,:])), norm_predictions[0,5])
        # Has the multiplication with the bins worked?
        assert np.array_equal(np.multiply(np.divide(raw_predictions[0, 5, :],np.sum(raw_predictions[0, 5, :])), bins),
                              values[0,5, :])

        # Calculate mean value and errors
        mean_values = np.sum(values, axis=-1)
        variances = np.sum(
                        np.multiply(
                            np.square(
                                np.subtract(bins[None,None,:], mean_values[:,:,None])),
                                norm_predictions), axis=-1)
        sd = np.sqrt(variances)


        # Assert summation is correct
        assert np.array_equal(np.sum(values[0, 5, :]), mean_values[0, 5])
        # Asser variance is correct
        assert np.array_equal(np.sum((mean_values[0, 5]-bins)**2*norm_predictions[0,5]), variances[0,5])

        line = "{i} {j} 0 {m:.8f} {c:.8f} {sd:.8f}\n"
        data = []

        for i in range(nres-1):
            for j in range(i+1, nres):
                if info_key == "dist":
                    if (abs(i - j) > min_sep) and (contact_probabilities[i][j] > pthres):  # and (fasta_seq[i] != 'G') and (fasta_seq[j] != 'G'):
                        data.append((i+1, j+1, mean_values[i,j], contact_probabilities[i,j], sd[i,j]))
                else:
                    if (abs(i - j) > min_sep) and (contact_probabilities[i][j] > pthres):  # and (fasta_seq[i] != 'G') and (fasta_seq[j] != 'G'):
                        data.append((i+1, j+1, mean_values[i,j], contact_probabilities[i,j], sd[i,j]))
                    if info_key == "theta":  # Theta is assymmetric, add the other side of the diagonal if applicable
                        if (abs(j - i) > min_sep) and (contact_probabilities[j][i] > pthres):  # and (fasta_seq[i] != 'G') and (fasta_seq[j] != 'G'):
                            data.append((j+1, i+1, mean_values[i,j], contact_probabilities[j,i], sd[j,i]))

        # Sort on the probability of contact
        data.sort(key=lambda x: x[3], reverse=True)

        content = [header]

        # Format each line according to CASP rr format
        for i, j, m, c, sd in data:
            content.append(line.format(i=i, j=j, m=m, c=c, sd=sd))
        content.append("END\n")
        if len(out_base_path) > 0:
            out_file = os.path.join(out_base_path, os.path.basename(input_file[:-4]) + bin_values.ending)
        else:
            out_file = input_file[:-4] + bin_values.ending
        with open(out_file, 'w') as contacts_handle:
            contacts_handle.write(''.join(content))
        # return ''.join(content)


def _virtual_cb_vector(residue):
    # get atom coordinates as vectors
    n = residue['N'].get_vector()
    c = residue['C'].get_vector()
    ca = residue['CA'].get_vector()
    # center at origin
    n = n - ca
    c = c - ca
    # find rotation matrix that rotates n -120 degrees along the ca-c vector
    rot = rotaxis(-pi*120.0/180.0, c)
    # apply rotation to ca-n vector
    cb_at_origin = n.left_multiply(rot)
    # put on top of ca atom
    cb = cb_at_origin + ca
    return cb


def pdb_to_npz(npz_name, pdb_file=False, mmCIF_file=False, std=1):
    """
    Convert a pdb/mcif to trRosetta distances/angles
    """

    if pdb_file:
        from Bio.PDB.PDBParser import PDBParser
        bio_parser = PDBParser(PERMISSIVE=1)
        structure_file = pdb_file
        structure_id = pdb_file.name[:-4]
    elif mmCIF_file:
        from Bio.PDB.MMCIFParser import MMCIFParser
        bio_parser = MMCIFParser()
        structure_file = mmCIF_file
        structure_id = mmCIF_file.name[:-4]
    else:
        print("No file given: one pdb or one mmCIF file has to be definied")
        sys.exit()

    # Load structure
    structure = bio_parser.get_structure(structure_id, structure_file)

    # Get residues and length of protein
    residues = []
    for chain in structure[0]:
        for residue1 in structure[0][chain.id]:
            if not is_aa(residue1):
                continue
            residues.append(residue1.get_resname())
    plen = len(residues)

    # Setup bins and step for the final matrix
    DIST_STEP = 0.5
    OMEGA_STEP = 15
    THETA_STEP = 15
    PHI_STEP = 15

    z_per_bin = std/DIST_STEP
    z_step = z_per_bin/2
    angle_z_step = 1

    minvalue = 0.01
    dist_wanted_bins = (20 - 2)/DIST_STEP
    omega_wanted_bins = 360/OMEGA_STEP
    theta_wanted_bins = 360/THETA_STEP
    phi_wanted_bins = 180/PHI_STEP
    cumm_cutoff = 0.899
    # cb_lst = get_cb_coordinates(open(args.pdb_file, 'r'), "A")
    # contact_mat = get_cb_contacts(len(residues))
    dist_mat = np.full((plen, plen, 37), minvalue/36)
    omega_mat = np.full((plen, plen, 25), minvalue/24)
    theta_mat = np.full((plen, plen, 25), minvalue/24)
    phi_mat = np.full((plen, plen, 13), minvalue/12)

    dist_bins = [i for i in np.arange(2, 2 + (dist_wanted_bins)*0.5, 0.5)]
    omega_bins = [i for i in np.arange(-180, -180 +
                  (omega_wanted_bins)*OMEGA_STEP, OMEGA_STEP)]
    theta_bins = [i for i in np.arange(-180, -180 +
                  (theta_wanted_bins)*THETA_STEP, THETA_STEP)]
    phi_bins = [i for i in np.arange(0, (phi_wanted_bins)*PHI_STEP, PHI_STEP)]
    dist_num_bins = len(dist_bins)
    omega_num_bins = len(omega_bins)
    theta_num_bins = len(theta_bins)
    phi_num_bins = len(phi_bins)

    # Iterate over all residues and calculate distances
    i = 0
    j = 0
    for chain in structure[0]:
        for residue1 in structure[0][chain.id]:
        # Only use real atoms, not HET or water
            if not is_aa(residue1):
                continue
            
            # If the residue lacks CB (Glycine etc), create a virtual
            if residue1.has_id('CB'):
                c1B = residue1['CB'].get_vector()
            else:
                c1B = _virtual_cb_vector(residue1)
            
            j = 0
            for chain in structure[0]:
                for residue2 in structure[0][chain.id]:
                    symm = False
                    if not is_aa(residue2):
                        continue
                    # print(i,j)
                    if i == j:
                        dist_mat[i, j, 0] = (1-minvalue)
                        omega_mat[i, j, 0] = (1-minvalue)
                        theta_mat[i, j, 0] = (1-minvalue)
                        phi_mat[i, j, 0] = (1-minvalue)
                        j += 1
                        continue
                    
                    if i > j:
                        dist_mat[i, j] = dist_mat[j, i]
                        omega_mat[i, j] = omega_mat[j, i]
                        symm = True
                    # If the residue lacks CB (Glycine etc), create a virtual
                    if residue2.has_id('CB'):
                        c2B = residue2['CB'].get_vector()
                    else:
                        c2B = _virtual_cb_vector(residue2)
                    ###############################################
                    dist = (c2B-c1B).norm()
                    
                    if dist > 20:
                        dist_mat[i, j, 0] = (1-minvalue)
                        omega_mat[i, j, 0] = (1-minvalue)
                        theta_mat[i, j, 0] = (1-minvalue)
                        phi_mat[i, j, 0] = (1-minvalue)
                    else:
                        # Dist and omega are symmetrical and have already been copied
                        if not symm:
                            ix = np.digitize(dist, dist_bins)
                            cum_prob = 0
                            b_step = 0
                            while cum_prob < cumm_cutoff:
                                bin_prob = st.norm.cdf(b_step*-z_step) -\
                                    st.norm.cdf(-z_step*(1+b_step))
                                dist_mat[i, j, np.min([ix+b_step, dist_num_bins])]\
                                    += bin_prob
                                dist_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                                cum_prob += bin_prob*2
                                b_step += 1
                        ###############################################
                        # # Omega
                        c1A = residue1['CA'].get_vector()
                        c2A = residue2['CA'].get_vector()
                    
                        if not symm:
                            raw_omega = calc_dihedral(c1A, c1B, c2B, c2A)
                            omega = (raw_omega*180)/pi
                    
                            ix = np.digitize(omega, omega_bins)
                            cum_prob = 0
                            b_step = 0
                            while cum_prob < cumm_cutoff:
                                bin_prob = st.norm.cdf(b_step*-angle_z_step) -\
                                    st.norm.cdf(-angle_z_step*(1+b_step))
                                omega_mat[i, j, np.min([ix+b_step, omega_num_bins])]\
                                    += bin_prob
                                omega_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                                cum_prob += bin_prob*2
                                b_step += 1
                    
                        ###############################################
                        # # Theta
                        N1 = residue1['N'].get_vector()
                    
                        raw_theta = calc_dihedral(N1, c1A, c1B, c2B)
                        theta = (raw_theta*180)/pi
                    
                        ix = np.digitize(theta, theta_bins)
                        cum_prob = 0
                        b_step = 0
                        while cum_prob < cumm_cutoff:
                            bin_prob = st.norm.cdf(b_step*-angle_z_step) -\
                                st.norm.cdf(-angle_z_step*(1+b_step))
                            theta_mat[i, j, np.min([ix+b_step, theta_num_bins])]\
                                += bin_prob
                            theta_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                            cum_prob += bin_prob*2
                            b_step += 1
                    
                        ###############################################
                        # # Phi
                    
                        raw_phi = calc_angle(c1A, c1B, c2B)
                        phi = (raw_phi*180)/pi
                    
                        ix = np.digitize(phi, phi_bins)
                        cum_prob = 0
                        b_step = 0
                        while cum_prob < cumm_cutoff:
                            bin_prob = st.norm.cdf(b_step*-angle_z_step) -\
                                st.norm.cdf(-angle_z_step*(1+b_step))
                            phi_mat[i, j, np.min([ix+b_step, phi_num_bins])]\
                                += bin_prob
                            phi_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                            cum_prob += bin_prob*2
                            b_step += 1
                    j += 1
            i += 1

    np.savez_compressed(npz_name,
                        dist=dist_mat,
                        omega=omega_mat,
                        theta=theta_mat,
                        phi=phi_mat)
