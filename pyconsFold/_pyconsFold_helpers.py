#!/usr/bin/env python3
import os
import stat
import sys
import argparse
import time
import glob
import re
import shutil
import subprocess
import tempfile
try:
    from .utils import npz_to_casp
    extra_compatible = True
except ModuleNotFoundError:
    extra_compatible = False


AA3TO1 = {"ALA": "A", "ASN": "N", "CYS": "C", "GLN": "Q", "HIS": "H",
          "LEU": "L", "MET": "M", "PRO": "P", "THR": "T", "TYR": "Y",
          "ARG": "R", "ASP": "D", "GLU": "E", "GLY": "G", "ILE": "I",
          "LYS": "K", "PHE": "F", "SER": "S", "TRP": "W", "VAL": "V"}
AA1TO3 = {v: k for k, v in AA3TO1.items()}

CBATOM = {"A": "cb", "N": "cb", "C": "cb", "Q": "cb", "H": "cb", "L": "cb",
          "M": "cb", "P": "cb", "T": "cb", "Y": "cb", "R": "cb", "D": "cb",
          "E": "cb", "G": "ca", "I": "cb", "K": "cb", "F": "cb", "S": "cb",
          "W": "cb", "V": "cb"}


def check_programs(program_dssp, program_pcons, program_tmscore, cns_suite, cns_executable):
        ### Start by checking so CNS is installed ###
        if not os.path.isdir(cns_suite):
            print("ERROR! Can not find CNS suite folder at " +
                  "{}".format(cns_suite))
            print("It does not look like you have CNS installed and in your path.")
            print("Install cns (http://cns-online.org/v1.3/) and follow the setup instructions (https://github.com/johnlamb/pyconsfold)")
            sys.exit()
        if not os.path.isfile(cns_suite + "/cns_solve_env.sh"):
            print("ERROR! cns_solve_env.sh not found inside " +
                  "{}".format(cns_suite))
            print("Check CNS installation!")
            sys.exit()
        if not os.path.isfile(cns_executable):
            print("ERROR! CNS executable not found at " +
                  "{}".format(cns_executable))
            print("Check CNS installation!")
            sys.exit()
        if not os.path.isfile(cns_suite + "/libraries/toppar/protein.param"):
            print("ERROR! protein.param not found inside " +
                  "{}/libraries/toppar/!".format(cns_executable))
            print("Check CNS installation!")
            sys.exit()

        ### Check the different utils ###
        if shutil.which("dssp"):
            program_dssp = shutil.which("dssp")
        elif not os.path.isfile(program_dssp):
            print("ERROR! Can not find dssp program at " +
                  "{}".format(program_dssp))
            sys.exit()

        if shutil.which("pcons"):
            program_dssp = shutil.which("pcons")
        elif not os.path.isfile(program_pcons):
            print("ERROR! Can not find pcons program at " +
                  "{}".format(program_pcons))
            sys.exit()

        if shutil.which("TMscore"):
            program_tmscore = shutil.which("TMscore")
        elif not os.path.isfile(program_tmscore):
            print("ERROR! You can only compare the models against a native structure if you have TMscore installed")
            print("ERROR! Can not find TMscore program at " +
                  "{}".format(program_TMscore))
            print("Please install TMscore to use")
            sys.exit()
        return program_dssp, program_pcons, program_tmscore, cns_suite, cns_executable


def seq_fasta(fasta_file):
    seq = ""
    with open(fasta_file) as fasta_handle:
        for line in fasta_handle:
            if line.startswith(">"):
                continue
            else:
                seq += line.strip()
    return seq


def seq_rr(rr_file, dir_out):
    seq_rr = ""
    with open(rr_file) as rr_handle:
        for line in rr_handle:
            if line.startswith("PFRMAT"):
                continue
            elif line.startswith("TARGET"):
                continue
            elif line.startswith("AUTHOR"):
                continue
            elif line.startswith("SCORE"):
                continue
            elif line.startswith("REMARK"):
                continue
            elif line.startswith("METHOD"):
                continue
            elif line.startswith("MODEL"):
                continue
            elif line.startswith("PARENT"):
                continue
            elif line.startswith("TER"):
                break
            elif line.startswith("END"):
                break
            elif re.match("^[0-9 -\.]+$", line):
                break
            else:
                seq_rr += line.strip()
    if len(seq_rr) == 0:
        print("No sequence header in {}".format(rr_file))
        clean_output_dir(dir_out)
        sys.exit()
    return seq_rr


def print2file(write_file, text, newline='\n'):
    if os.path.isfile(write_file):
        with open(write_file, 'a+') as write_handle:
            write_handle.write('\n' + text)
    else:
        with open(write_file, 'w') as write_handle:
            write_handle.write(text)


def flatten_fasta(fasta_file):
    seq = seq_fasta(fasta_file)
    header = open(fasta_file).readline().strip()
    os.remove(fasta_file)
    print2file(fasta_file, header + '\n' + seq)


def fasta2residues(fasta_file):
    seq = seq_fasta(fasta_file)
    residues = {i+1: res for i, res in enumerate(seq)}
    return residues


def rr2contacts(rr_file, seq_sep):
    contacts = {}
    with open(rr_file) as rr_handle:
        for line in rr_handle:
            if not re.match("^[0-9]", line):
                continue
            else:
                c = [float(x) for x in line.split()]
                if abs(c[1]-c[0]) < seq_sep:
                    continue
                else:
                    # Compatability for both with and without distance errors
                    if len(c)>5:
                        cont = [c[3], c[4], c[5]]
                    else:
                        cont = [c[3], c[4]]
                    contacts[" ".join([str(int(c[0])),
                                       str(int(c[1]))])] = cont
    return contacts


def angle2restraints(angle_file, name, seq_sep):
    contacts = {}
    with open(angle_file) as angle_handle:
        for line in angle_handle:
            if not re.match("^[0-9]", line):
                continue
            else:
                c = [float(x) for x in line.split()]
                if abs(c[1]-c[0]) < seq_sep:
                    continue
                else:
                    # contacts.append([int(x) for x in c[:2]])
                    # Compatability for both with and without distance errors
                    cont = [c[3], c[4], c[5]]
                    contacts[" ".join([str(int(c[0])),
                                       str(int(c[1])), name])] = cont
    return contacts


def chk_errors_seq(seq, dir_out):
    for res in seq:
        if res not in AA1TO3:
            print("Undefined amino acid {} in {}".format(res, seq))
            clean_output_dir(dir_out)
            sys.exit()


def clean_output_dir(dir_out):
    # for fn in glob.glob(dir_out + "/stage1/*_model*.pdb"):
    #     shutil.copy(fn, dir_out)
    shutil.rmtree(dir_out + "/input")
    if os.path.isdir(dir_out + "/stage1"):
        shutil.rmtree(dir_out + "/stage1")
    shutil.rmtree(dir_out + "/stage1")
    if os.path.isdir(dir_out + "/stage2"):
        shutil.rmtree(dir_out + "/stage2")


# pair not implemented
def process_arguments(fasta, contacts, dir_out, ss, rrtype, selectrr, fasta2, omega, theta, mcount,
                      lbd, contwt, sswt, rep2, rr_pthres, rr_sep, use_angles, tmscore_pdb_file, debug):
    if tmscore_pdb_file:
        if not os.path.isfile(tmscore_pdb_file):
            print("ERROR! Native PDB file {} does not exist!".format(tmscore_pdb_file))
            sys.exit()

    npz=False
    rr=False
    if contacts.lower().endswith(".npz"):
        npz=True
        ## If we use a npz file, default rr_pthres to 0.55 as we have n^2 contacts(and because benchmark shows it)
        if rr_pthres == 0:
            rr_pthres = 0.55
    elif contacts.lower().endswith(".rr"):
        rr=True
    ### Check that all files exists ###
    if not os.path.isfile(fasta):
        print("ERROR! Fasta file {} does not exist!".format(fasta))
        sys.exit()
    if npz:
        if not os.path.isfile(contacts):
            print("ERROR! Contact file {} does not exist!".format(contacts))
            sys.exit()
    elif rr:
        if not os.path.isfile(contacts):
            print("ERROR! Contact file {} does not exist!".format(contacts))
            sys.exit()
    else:
        print("ERROR! At least one contacts file needs to be supplied with either .rr (CASP) or .npz (trRosetta) file-ending")
        sys.exit()

    if fasta2:
        if not os.path.isfile(fasta2):
            print("ERROR! Fasta file #2 {} does not exist!".format(fasta2))
            sys.exit()
    if ss:
        if not os.path.isfile(ss):
            print("ERROR! Secondary structure file {} " +
                  "does not exist!".format(fasta))
            sys.exit()

    ### Length and base name of fasta seq ###
    L = len(seq_fasta(fasta))
    mini = 15 * L
    f_id = os.path.splitext(os.path.basename(fasta))[0]
    f2_id = os.path.splitext(os.path.basename(fasta2))[0]
    if f2_id == f_id:
        f2_id += '_2'

    ### How many contacts to use? ###
    selectrr = selectrr.replace("L", "")
    selectrr = 10000000 if selectrr == "all" else int(float(selectrr) * L + 0.5)
    if not (selectrr > 0 and selectrr <= 10000000):
        print("ERROR! Selectrr option does not " +
              "look right: {}".format(selectrr))
        sys.exit()

    ### mcount, how many models to generate ###
    if not (mcount >= 5 and mcount <= 50):
        print("ERROR! Model count must be between 5 and 50!")
        sys.exit()

    ### Lambda value ###
    if not (lbd >= 0.1 and lbd <= 10.0):
        print("ERROR! Lambda must be between 0.1 and 10.0!")
        sys.exit()

    ### Contact restraint weight, higher means more weight to contact restraints ###
    if not (contwt >= 0.1 and contwt <= 10000):
        print("ERROR! Contact restraint weights must be " +
              "between 0.1 and 10000!")
        sys.exit()

    ### Secondary structure weight, more means more weight to secondary structure restraints ###
    if not (sswt >= 0.1 and sswt <= 100):
        print("ERROR! SS restraint weights must be between 0.1 and 100!")
        sys.exit()

    ### Rep2 ###
    if not (rep2 >= 0.6 and rep2 <= 1.5):
        print("ERROR! Rep2 must be between 0.6 and 1.5!")
        sys.exit()


    ### Prepare working directories, remove if exists ###
    dir_out = os.path.abspath(os.path.join(os.getcwd(), dir_out))
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    else:
        shutil.rmtree(dir_out)
        os.makedirs(dir_out)
    if os.path.isdir(dir_out + "/input"):
        shutil.rmtree(dir_out + "/input")
    if os.path.isdir(dir_out + "/stage1"):
        shutil.rmtree(dir_out + "/stage1")
    if os.path.isdir(dir_out + "/stage2"):
        shutil.rmtree(dir_out + "/stage2")
    os.makedirs(dir_out + "/input")

    ## Create filenames to standardize ###
    fasta_file = f_id + ".fasta"
    rr_file = f_id + ".rr"

    ### Optional files ###
    if npz:
        if not extra_compatible:
            print("ERROR: Using npz-files requires additional dependencies (BioPython, scypi and numpy)")
            print("Install the using: `pip3 install -U pyconsFold[all]`")
            clean_output_dir(dir_out)
            sys.exit()
        npz_dir = tempfile.mkdtemp()
        if fasta2:
            npz_to_casp(contacts, fasta_file=fasta, fasta2_file=fasta2, out_base_path=npz_dir)
        else:
            npz_to_casp(contacts, fasta_file=fasta, out_base_path=npz_dir)
        contacts_file = os.path.join(npz_dir, os.path.basename(contacts[:-4]) + '.rr')
        if use_angles:
            omega = os.path.join(npz_dir, os.path.basename(contacts[:-4]) + '.omega')
            theta = os.path.join(npz_dir, os.path.basename(contacts[:-4]) + '.theta')
            phi = os.path.join(npz_dir, os.path.basename(contacts[:-4]) + '.phi')
    else:
        contacts_file = contacts
    if fasta2:
        fasta2_file = f2_id + ".fasta"
    else:
        fasta2_file = ''
    if ss:
        ss_file = f_id + ".ss"
    else:
        ss_file = ''
    if omega:
        omega_file = f_id + ".omega"
    else:
        omega_file = ''
    if theta:
        theta_file = f_id + ".theta"
    else:
        theta_file = ''

    ### Copy files into the input directory ###
    shutil.copy(fasta, dir_out + "/input/" + fasta_file)
    shutil.copy(contacts_file, dir_out + "/input/" + rr_file)
    if fasta2:
        shutil.copy(fasta2, dir_out + "/input/" + fasta2_file)
    if ss:
        shutil.copy(ss, dir_out + "/input/" + ss_file)
    if omega:
        shutil.copy(omega, dir_out + "/input/" + omega_file)
    if theta:
        shutil.copy(theta, dir_out + "/input/" + theta_file)
    if npz:
        shutil.rmtree(npz_dir)

    ### Move into the input directory ###
    # base_dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.getcwd()
    os.chdir(dir_out + "/input")

    ### Flatten fasta, make sure it has the sequence on one row ###
    flatten_fasta(fasta_file)

    ### Check integrity of fasta file ###
    residues = fasta2residues(fasta_file)
    seq = seq_fasta(fasta_file)
    plen = len(seq)
    p2len = 0
    chk_errors_seq(seq, dir_out)

    ### If we use a second fasta (for docking), check for error and concatenate sequences ###
    if fasta2:
        flatten_fasta(fasta2_file)
        seq2 = seq_fasta(fasta2_file)
        p2len = len(seq2)
        chk_errors_seq(seq2, dir_out)

        ### Update both seq and residues to use both files
        fasta2_res = fasta2residues(fasta2_file) 
        inter_dict = {}
        for key, value in fasta2_res.items():
            inter_dict[key+plen] = value
        residues.update(inter_dict)
        seq = seq + seq2

    ### Check so that the sequence from the contact file is the same as the sequence from the fasta file ###
    ### (Must be done after a potential second file has been loaded ###
    rr_seq = seq_rr(rr_file, dir_out)
    if not seq == rr_seq:
        print("ERROR! Fasta and rr sequence do not match!" +
              "\nFasta\t: {} \nRR\t: {}".format(seq, rr_seq))
        clean_output_dir(dir_out)
        sys.exit()

    ### Check SS file for integrity and length against fasta
    if ss:
        ss_seq = seq_fasta(ss_file)
        if len(seq) != len(ss_seq):
            print("ERROR! Fasta and ss sequence length do not match!" +
                  "\nFasta\t: {} \nRR\t:{}".format(seq, ss_seq))
            clean_output_dir(dir_out)
            sys.exit()

        for s in ss_seq:
            if s not in "HCE":
                print("ERROR undefined secondary structure unit {}".format(s))
                clean_output_dir(dir_out)
                sys.exit()

    ### Read in contacts and sort them on probability (5th column) ### 
    rr_lines = []
    with open(rr_file) as rr_handle:
        contact_between_proteins = False
        for line in rr_handle:
            if re.match('^[A-Za-z]', line):
                continue
            if len(line.strip()) == 0:
                continue
            x, y = (int(x) for x in line.strip().split()[:2])
            contact_prob = float(line.strip().split()[4])
            ### Only use contacts according to separation ###
            if (y-x)>rr_sep and contact_prob > rr_pthres:
                rr_lines.append(line.strip().split())
                if fasta2:
                    #  If we do docking, make sure we at least have one contact between the proteins
                    if ((x < (plen+1)) and (y > plen)) or ((y < (plen+1)) and (x > plen)):
                        contact_between_proteins = True

    if fasta2 and not contact_between_proteins:
        # We are doing docking but there are no contacts between proteins, adding an artificial one
        # Make a contact between the middle of each protein with distance (len(A)/2+len(B)/2)*3.5Å
        print("WARNING, trying to do docking without any predicted contacts between the proteins.")
        print("Adding a weak and distant artificial contact")
        print("This increases the modelling time significantly")
        mid_first = plen//2
        mid_second = p2len//2
        a_dist = (mid_first + mid_second)*3.5
        print("From {} to {} with distance {}Å and error {}Å".format(mid_first,plen+mid_second,a_dist, a_dist/4))
        rr_lines.append([str(mid_first), str(plen+mid_second),'0',str(a_dist),str(rr_pthres+0.01),str(a_dist/4)])
    ### Sort and only use contacts above set threshold ###
    rr_lines.sort(key=lambda x: float(x[4]), reverse=True)
    rr_scores = '\n'.join([' '.join(x) for x in rr_lines if float(x[4]) > rr_pthres])

    ### Remove old contacts file and replace with new sorted and filtered ###
    os.remove(rr_file)
    print2file(rr_file, seq + '\n' + rr_scores + '\n')

    contacts = rr2contacts(rr_file, 1)

    ### Make sure contacts and residues are synced ###
    for key in contacts.keys():
        a, b = key.split()
        if not residues[int(a)]:
            print("ERROR! Residue {} is out of sequence!".format(a))
            clean_output_dir(dir_out)
            sys.exit()
        if not residues[int(b)]:
            print("ERROR! Residue {} is out of sequence!".format(b))
            clean_output_dir(dir_out)
            sys.exit()

    ### Change back to base directory ###
    os.chdir(base_dir)

    ### If debug is on, print all variables ###
    if debug:
        print("dir_out      {}".format(dir_out))
        print("fasta_file   {}".format(fasta_file))
        print("rr_file      {}".format(rr_file))
        print("ss_file      {}".format(ss_file))
        print("omega_file   {}".format(omega_file))
        print("theta_file   {}".format(theta_file))
        print("selectrr     {}".format(selectrr))
        print("rrtype       {}".format(rrtype))
        print("Lambda       {}".format(lbd))
        print("id           {}".format(f_id))
        print("L            {}".format(L))
        print("sequence     {}".format(seq))
        if ss:
            print("ss_seq       {}".format(ss_seq))

    return fasta_file, fasta2_file, rr_file, ss_file, omega_file, theta_file, residues,\
        f_id, selectrr, mini


def write_cns_seq(fasta_file, file_cns_seq, chunk=64):
    with open(fasta_file) as fasta_handle:
        # Strip the header
        fasta_handle.readline().strip()
        seq = fasta_handle.readline().strip()
    three_letter_seq = ' '.join([AA1TO3[c] for c in seq])
    if os.path.isfile(file_cns_seq):
        os.remove(file_cns_seq)
    width_limited_seq = '\n'.join([three_letter_seq[i: chunk + i]
                                   for i in range(0,
                                                  len(three_letter_seq),
                                                  chunk)])
    print2file(file_cns_seq, width_limited_seq)


def build_extended(fasta_file, cns_suite, cns_executable, fasta_file2=""):
    write_cns_seq(fasta_file, "input.seq")
    if fasta_file2:
        input_len = len(seq_fasta(fasta_file))
        write_cns_seq(fasta_file2, "input2.seq")
        write_cns_gseq("input2.seq", input_len + 1)
    else:
        write_cns_gseq("", 0)

    # Not needing to generate gseq.inp file, already done as template
    # Not needing to generate extn.inp file, already done as template
    job_file = "#!/bin/bash\n"
    job_file += "source {}/cns_solve_env.sh\n".format(cns_suite)
    job_file += "export KMP_AFFINITY=none\n"
    job_file += "{} < gseq.inp > gseq.log\n".format(cns_executable)
    job_file += "{} < extn.inp > extn.log".format(cns_executable)
    print2file("job.sh", job_file)
    st = os.stat("job.sh")
    os.chmod("job.sh", st.st_mode | stat.S_IEXEC)
    # subprocess.call("chmod +x job.sh", shell=True)
    # subprocess.call("./job.sh &> job.log", shell=True)
    ######### No shell, does not work on Keb #############
    subprocess.call("./job.sh")
    if not os.path.isfile("extended.pdb"):
        print("FAILED! extended.pdb not found")
        sys.exit()
    os.remove("gseq.log")
    os.remove("extn.log")


def parse_pdb_row(row, param):
    switcher = {
        "anum": [6, 5],
        "aname": [12, 4],
        "altloc": [16, 1],
        "rname": [17, 3],
        "rnum": [22, 5],
        "insertion": [26, 1],
        "chain": [21, 1],
        "x": [30, 8],
        "y": [38, 8],
        "z": [46, 8]
        }
    offset, length = switcher[param]
    # print(offset, length)
    result = row[offset:offset + length].strip()
    return result


def pdb2rnum_rname(chain):
    rnum_rname_dict = {}
    with open(chain) as chain_handle:
        for line in chain_handle:
            if not line.startswith("ATOM"):
                continue
            else:
                rnum_rname_dict[parse_pdb_row(line, "rnum")] =\
                        parse_pdb_row(line, "rname")
    return rnum_rname_dict


def xyz_pdb(chain, atom_selection):
    atom_selection = atom_selection.upper()
    xyz_dict = {}
    if not os.path.isfile(chain):
        print("ERROR! File {} does not exists!".format(chain))
        sys.exit()
    if atom_selection not in ["CA", "CB", "ALL"]:
        print("ERROR! Selection must be ca, cb or all!")
        sys.exit()
    with open(chain) as chain_handle:
        for line in chain_handle:
            if not line.startswith("ATOM"):
                continue
            else:
                xyz_dict[parse_pdb_row(line, "rnum") + " " +
                         parse_pdb_row(line, "aname")] =\
                           " ".join([parse_pdb_row(line, "x"),
                                     parse_pdb_row(line, "y"),
                                     parse_pdb_row(line, "z")])

    if not xyz_dict:
        print("ERROR! xyz_pdb is empty")
        sys.exit()

    if atom_selection == "ALL":
        return xyz_dict

    rnum_rname = pdb2rnum_rname(chain)
    selected_xyz = {}
    for key in sorted(xyz_dict.keys()):
        C = key.split()
        this_atom = "CA" if (atom_selection == "CB" and
                             rnum_rname[C[0]] == "GLY") else atom_selection
        if C[1] != this_atom:
            continue
        else:
            selected_xyz[C[0]] = xyz_dict[key]

    return selected_xyz


def rr2r1a1r2a2(rr_file, rrtype, dir_out):
    contacts = rr2contacts(rr_file, 1)
    seq = seq_rr(rr_file, dir_out)
    r1a1r2a2 = {}
    for key in contacts.keys():
        a, b = [int(x) for x in key.split()]
        r1 = seq[a-1]
        r2 = seq[b-1]
        ca1 = ca2 = rrtype
        if rrtype.upper() == "CB":
            ca1 = CBATOM[r1]
            ca2 = CBATOM[r2]
        r1a1r2a2[" ".join([str(a), ca1, str(b), ca2])] = contacts[key]
    return r1a1r2a2


def rr2tbl(rr_file, tbl_file, rrtype, dir_out, dist):
    r1a1r2a2 = rr2r1a1r2a2(rr_file, rrtype, dir_out)
    if os.path.isfile(tbl_file):
        os.remove(tbl_file)
    to_print = ["nrestraints={}".format(len(r1a1r2a2) + 1)]
    for key, value in sorted(r1a1r2a2.items(), key=lambda i: i[1][1], reverse=True):
        # print(key, "=>", r1a1r2a2[key])
        C = key.split()
        if dist:
            if len(value) == 2:
                print("Warning! The contacts file lack a predicted error column,")
                print("using default -0.10Å and dist - 3.60Å as errors")
                negdev = 0.10
                posdev = value[0] - 3.60
            else:
                negdev = posdev = value[2]
        else:
            negdev = 0.10
            # Assuming less than 8Å is a contact
            posdev = 8 - 3.60
        if dist:
            distance = value[0]
        else:
            distance = 3.6
        to_print.append("assign (resid {:>3} and ".format(C[0]) +
                        "name {:>2}) ".format(C[1]) +
                        "(resid {:>3} and name {:>2}) ".format(C[2], C[3]) +
                        "{:.2f} {:.2f} ".format(distance, negdev) +
                        "{:.2f}".format(posdev))
    print2file(tbl_file, '\n'.join(to_print) + '\n')


def residue_atom_coordinates(xyz_dict):
    residues = {}
    for key, value in xyz_dict.items():
        res, atm = key.split()
        res = int(res)
        if res in residues:
            residues[res].append([atm, value])
        else:
            residues[res]=[[atm, value]]
    return residues


def contact_restraints(stage, selectrr, rrtype, dir_out, dist, f_id, rr_file=None, debug=False):
    if stage == "stage1":
        if not rr_file:
            return
        xL = selectrr
        # print(rr_file)
        rr_data = []
        with open(rr_file) as rr_handle:
            i = 0
            for line in rr_handle:
                rr_data.append(line)
                if not re.match("^[0-9]", line):
                    continue
                i += 1
                if i >= xL:
                    break
        os.remove(rr_file)
        print2file(rr_file, ''.join(rr_data))

    elif stage == "stage2":
        stg1_rr = os.path.join(dir_out, "stage1", f_id + ".rr")
        stg1_model = os.path.join(dir_out, "stage1", f_id + "_model1.pdb")
        xyz_dict = xyz_pdb(stg1_model, "ALL")
        res_atom_coords = residue_atom_coordinates(xyz_dict)

        rr_header = []
        rr_data = []
        rr_accepted = []
        rr_ignored = []
        with open(stg1_rr) as rr_handle:
            for line in rr_handle:
                if re.match("^[0-9]", line):
                    c = [float(x) for x in line.split()]
                    r1, r2 = c[:2]
                    pred_dist = c[3]
                    min_dist = 10000
                    for left in res_atom_coords[r1]:
                        for right in res_atom_coords[r2]:
                            measured_dist = calc_dist(left[1],right[1])
                            if measured_dist < min_dist:
                                min_dist = measured_dist
                    if min_dist < pred_dist:
                        rr_data.append(line)
                        rr_accepted.append("Accepted " + line.strip() + " min_dist={}".format(min_dist))
                    else:
                        rr_ignored.append("Ignored " + line.strip() + " min_dist={}".format(min_dist))

                else:
                    rr_header.append(line)
        if len(rr_data) == 0:
            print("ERROR!! No contacts accepted for stage2")
            sys.exit()
        rr_header.extend(rr_data)
        rr_accepted.extend(rr_ignored)
        print2file(rr_file, ''.join(rr_header))
        print2file(f_id + "_contact_filter.log", '\n'.join(rr_accepted))
        if debug:
            print("Stage 2 contact filtering.. removed {} contacts".format(len(rr_ignored)))

    rr2tbl(rr_file, "contact.tbl", rrtype, dir_out, dist)

                

def angle_restraints(omega_file, use_omega, theta_file, use_theta, residues, seq_sep=1):
    dihedral_file = "dihedral.tbl"
    # print(theta_file)
    angle_contacts = {} 
    # print(use_omega, use_theta)
    if use_omega:
        omega_contacts = angle2restraints(omega_file, 'omega', seq_sep)
        angle_contacts.update(omega_contacts)
        # print(len(angle_contacts))
    if use_theta:
        theta_contacts = angle2restraints(theta_file, 'theta', seq_sep)
        angle_contacts.update(theta_contacts)
        # print(len(angle_contacts))
    dihedral_to_write = []
    n = 0
    # mixed_contacts = sorted(omega_contacts.items(), key=lambda i: i[1][1], reverse=True)  # [:500]
    # mixed_contacts.extend(sorted(theta_contacts.items(), key=lambda i: i[1][1], reverse=True))  # [:500])
    # extend(sorted(theta_contacts.items(), key=lambda i: i[1][1], reverse=True)[:500])

    # Sort again, all angles mixed
    # print(len(angle_contacts))
    for key, value in sorted(angle_contacts.items(), key=lambda i: i[1][1], reverse=True):
        i, j, name = key.split()
        # # If one of the residues is Glycine, move along
        # if 'G' in (residues[int(i)] + residues[int(j)]):
        #     continue
        angle_mean, prob, angle_error = value

        if name == "omega":
            dihedral_to_write.append(("assign (resid {:>3} and name ca) " +
                                      "(resid {:>3} and name  cb) (resid" +
                                      " {:>3} and name cb) (resid {:>3}" +
                                      " and name ca) 1.0 {:>7} {:>7} 2")
                                      .format(i, i, j, j, angle_mean, angle_error))
        elif name == "theta":
            dihedral_to_write.append(("assign (resid {:>3} and name n) " +
                                      "(resid {:>3} and name  ca) (resid" +
                                      " {:>3} and name cb) (resid {:>3}" +
                                      " and name cb) 1.0 {:>7} {:>7} 2")
                                      .format(i, i, i, j, angle_mean, angle_error))
        n += 1
        # Only read in max 1500 dihedral angles, hard block in CNS
        if n > 1500:
            break

    # if res_sec:
    #     log_to_write = "write helix tbl restrains"
    #     dihedral_to_write = []
    #     hbnd_to_write = []
    #     for i in sorted(res_sec.keys()):
    #         PHI = res_dihe["H PHI"].split()
    #         PSI = res_dihe["H PSI"].split()
    #         if i-1 in residues.keys():
    #             dihedral_to_write.append(("assign (resid {:>3} and name c) " +
    #                                       "(resid {:>3} and name  n) (resid" +
    #                                       " {:>3} and name ca) (resid {:>3}" +
    #                                       " and name c) 5.0 {:>7} {:>7} 2 " +
    #                                       "!helix phi")
    #                                      .format(i-1, i, i, i, PHI[0], PHI[1]))
    #         if i+1 in residues.keys():
    #             dihedral_to_write.append(("assign (resid {:>3} and name n) " +
    #                                       "(resid {:>3} and name ca) (resid" +
    #                                       " {:>3} and name  c) (resid {:>3}" +
    #                                       " and name n) 5.0 {:>7} {:>7} 2 " +
    #                                       "!helix psi")
    #                                      .format(i, i, i, i+1, PSI[0], PSI[1]))
    print2file(dihedral_file, "\n".join(dihedral_to_write))


def print_helix_restraints(ss_file, residues, log_file, res_dihe,
                           res_hbnd, res_dist, ATOMTYPE, SHIFT):
    hbnd_file = "hbond.tbl"
    dihedral_file = "dihedral.tbl"
    ssnoe_file = "ssnoe.tbl"
    res_sec = fasta2residues(ss_file)
    res_sec = {key: value for key, value in res_sec.items() if value == "H"}

    if res_sec:
        log_to_write = "write helix tbl restrains"
        dihedral_to_write = []
        hbnd_to_write = []
        for i in sorted(res_sec.keys()):
            PHI = res_dihe["H PHI"].split()
            PSI = res_dihe["H PSI"].split()
            if i-1 in residues.keys():
                dihedral_to_write.append(("assign (resid {:>3} and name c) " +
                                          "(resid {:>3} and name  n) (resid" +
                                          " {:>3} and name ca) (resid {:>3}" +
                                          " and name c) 5.0 {:>7} {:>7} 2 " +
                                          "!helix phi")
                                         .format(i-1, i, i, i, PHI[0], PHI[1]))
            if i+1 in residues.keys():
                dihedral_to_write.append(("assign (resid {:>3} and name n) " +
                                          "(resid {:>3} and name ca) (resid" +
                                          " {:>3} and name  c) (resid {:>3}" +
                                          " and name n) 5.0 {:>7} {:>7} 2 " +
                                          "!helix psi")
                                         .format(i, i, i, i+1, PSI[0], PSI[1]))
        print2file(dihedral_file, "\n".join(dihedral_to_write))

        for i in sorted(res_sec.keys()):
            if i+1 not in res_sec.keys():
                continue
            elif i+2 not in res_sec.keys():
                continue
            elif i+3 not in res_sec.keys():
                continue
            elif i+4 not in res_sec.keys():
                continue
            HATOM = "H"
            HR = [float(x) for x in res_hbnd[HATOM].split()]
            # In case of Proline, use N atom instead of H
            HR[0] -= 1.0
            if residues[i+4] == "P":
                HR[0] += 1.0
                HATOM = "N"
            hbnd_to_write.append(("assign (resid {:>3} and name O) (resid " +
                                  "{:>3} and name {}) {:>4.2f} {:>4.1f} " +
                                  "{:>4.1f} !helix")
                                 .format(i, i+4, HATOM, HR[0], HR[1], HR[2]))
        if len(hbnd_to_write) > 0:
            print2file(hbnd_file, "\n".join(hbnd_to_write) + '\n')

        ssnoe_to_write = []
        for i in sorted(res_sec.keys()):
            if i+1 not in res_sec.keys():
                continue
            elif i+2 not in res_sec.keys():
                continue
            elif i+3 not in res_sec.keys():
                continue
            elif i+4 not in res_sec.keys():
                continue
            for A in sorted(ATOMTYPE.keys()):
                for SH in sorted(SHIFT.keys()):
                    AR = [float(x) for x in res_dist["H {}-{} O {}"
                                                     .format(A, A, SH)]
                          .split()]
                    if i+4+int(SH) not in res_sec.keys():
                        continue
                    ssnoe_to_write.append(("assign (resid {:>3} and name " +
                                           "{:>2}) (resid {:>3} and name " +
                                           "{:>2}) {:.2f} {:.2f} {:.2f} " +
                                           "!helix")
                                          .format(i, A, i+4+int(SH),
                                                  A, AR[0], AR[1], AR[2]))
        if len(ssnoe_to_write) > 0:
            print2file(ssnoe_file, "\n".join(ssnoe_to_write))

    else:
        log_to_write = "no helix predictions!"

    print2file(log_file, log_to_write)


def strand_count(ss_file):
    ss_seq = seq_fasta(ss_file)
    count = 0
    for i in range(1, len(ss_seq) - 1):
        if ss_seq[i] != 'E':
            continue
        elif ss_seq[i-1] == 'E' and ss_seq[i+1] != 'E':
            count += 1
    return count


def pairing2hbonds(pair_file, pair_hbonds):
    # TODO
    pass


def strand_and_sheet_tbl(ss_file, pair_file, pair_hbonds,
                         res_dihe, res_strnd_OO):
    res_ss = fasta2residues(ss_file)
    res_ssE = res_ss.copy()
    res_ssE = {key: value for key, value in res_ssE.items() if value == "E"}
    if pair_file:
        # TODO Implement together with pairing2hbonds
        pass

    # Identify the strands that are not used for pairing and generate
    # generic dihedral restraints for them
    dihedral_file = "dihedral.tbl"
    dihedral_to_write = []
    ssnoe_file = "ssnoe.tbl"
    ssnoe_to_write = []
    for i in sorted(res_ssE.keys()):
        strand_type = "unpaired E residue"
        if pair_file:
            # TODO Implement together with pairing2hbonds
            pass
        else:
            SPHI = res_dihe["U PHI"].split()
            SPSI = res_dihe["U PSI"].split()

        if i-1 in res_ss.keys() and res_ss[i-1] == "E":
            dihedral_to_write.append(("assign (resid {:>3} and name c) " +
                                      "(resid {:>3} and name  n) (resid " +
                                      "{:>3} and name ca) (resid {:>3} and " +
                                      "name c) 5.0 {:>7} {:>7} 2 !{} phi")
                                     .format(i-1, i, i, i, SPHI[0], SPHI[1],
                                             strand_type))
        if i+1 in res_ss.keys() and res_ss[i+1] == "E":
            dihedral_to_write.append(("assign (resid {:>3} and name n) " +
                                      "(resid {:>3} and name ca) (resid " +
                                      "{:>3} and name  c) (resid {:>3} and " +
                                      "name n) 5.0 {:>7} {:>7} 2 !{} psi")
                                     .format(i, i, i, i+1, SPSI[0], SPSI[1],
                                             strand_type))
    print2file(dihedral_file, "\n".join(dihedral_to_write) + '\n')
    for i in sorted(res_ssE.keys()):
        strand_type = "unpaired E residue"
        SD = []
        if pair_file:
            # TODO implement pair file check
            print("ERROR! Pair file not implemented yet")
            sys.exit()
        else:
            SD = [float(x) for x in res_strnd_OO["U"].split()]
        if i+1 not in res_ssE.keys():
            continue
        elif res_ssE[i+1] != "E":
            continue
        ssnoe_to_write.append(("assign (resid {:>3} and name {:>2}) (resid " +
                               "{:>3} and name {:>2}) {:.2f} {:.2f} {:.2f} " +
                               "!{}")
                              .format(i, "O", i+1, "O",
                                      SD[0], SD[1], SD[2], strand_type))
    print2file(ssnoe_file, "\n".join(ssnoe_to_write) + '\n')


def sec_restraints(stage, ss_file, res_dihe, res_hbnd, res_dist,
                   res_strnd_OO, residues, ATOMTYPE, SHIFT, debug):
    log_file = "ssrestraints.log"
    if os.path.isfile(log_file):
        os.remove(log_file)
    # helix restraints; start writing to the files
    # hbond.tbl, ssnoe.tbl and dihedral.tbl
    print_helix_restraints(ss_file, residues, log_file, res_dihe,
                           res_hbnd, res_dist, ATOMTYPE, SHIFT)
    strands = strand_count(ss_file)
    if strands == 0:
        print2file(log_file, "no strand restraints")
        return
    elif strands == 1:
        if debug:
            print("WARNING! Only 1 strand! Check ss-file")
        print2file(log_file, "WARNING! only 1 strand! please check your " +
                             "secondary structure file!")
        return

    if stage == "stage1":
        # if pair_file:
        #     ###################
        #     # NOT IMPLEMENTED #
        #     ###################
        #     pairing2hbonds(pair_file, "pairing_hbonds.txt")
        # else:
        strand_and_sheet_tbl(ss_file, None, None, res_dihe, res_strnd_OO)
    else:
        # TODO implement the next stage, using stage1 model
        # Now it only cares for helices, beta pairing is NOT implemented
        strand_and_sheet_tbl(ss_file, None, None, res_dihe, res_strnd_OO)


def update_nrestraints():
            ### Update the contact.tbl file with new number of restraints
            extra_restraints = 0
            for tbl in ["hbond.tbl", "dihedral.tbl", "ssnoe.tbl"]:
                if os.path.isfile(tbl):
                    with open(tbl,'r') as tbl_handle:
                        num = len(tbl_handle.read().split('\n'))
                        extra_restraints += num
            contact_to_write = ''
            with open("contact.tbl", 'r') as contact_handle:
                curr_restraints = int(contact_handle.readline().split('=')[1])
                contact_to_write += 'nrestraints=' + str(extra_restraints + curr_restraints) + '\n'
                for line in contact_handle:
                    contact_to_write += line
            os.remove("contact.tbl")
            print2file("contact.tbl", contact_to_write)


def count_lines(input_file):
    if not os.path.isfile(input_file):
        return 0
    lines = 0
    with open(input_file) as input_handle:
        for line in input_handle:
            w_line = line.strip()
            if len(w_line) < 1:
                continue
            lines += 1
    if input_file.endswith("contact.tbl"):
        lines -= 1
    return lines


def tbl2rows(tbl_file):
    noe = {}
    with open(tbl_file) as tbl_handle:
        for line in tbl_handle:
            w_line = line.strip()
            # This also handles the contact.tbl first line(setting nrestraints)
            if not w_line.startswith("assign"):
                continue
            w_line = re.sub("[\(\)]", " ", w_line)
            noe[w_line] = 1
    if not noe:
        print("{} seems empty!".format(tbl_file))
    return noe


def coverage_tbl(fasta_file, tbl, flag_dihe, fasta_file2=''):
    seq = seq_fasta(fasta_file)
    L = len(seq)
    if len(fasta_file2)>0:
        seq += seq_fasta(fasta_file2)
        L = len(seq)
    cov = re.sub("[A-Z]", "-", seq)
    noe = tbl2rows(tbl)

    for key in noe.keys():
        # assign (resid 123 and name CA) (resid  58 and name CA)
        # 3.60 0.10 3.40
        C = key.split()
        print(C)
        print(cov)
        r1 = int(C[2])
        r2 = int(C[7])
        # the case when "ca or cb" restraints are provided
        # assign ((resid 123 and name ca) or (resid 123 and name cb))
        # ((resid 58 and name ca) or (resid 58 and name cb)) 3.60 0.10 3.40
        if C[6] == "or":
            r2 = int(C[13])
        if flag_dihe:
            r2 = int(C[17])
        c1 = cov[r1-1]
        c2 = cov[r2-1]
        if c1 == "-":
            c1 = 1
        elif c1 == "*":
            c1 = "*"
        else:
            c1 = int(c1) + 1
            c1 = "*" if c1 > 9 else c1

        if c2 == "-":
            c2 = 1
        elif c2 == "*":
            c2 = "*"
        else:
            c2 = int(c2) + 1
            c2 = "*" if c2 > 9 else c2

        cov = cov[:r1-1] + str(c1) + cov[r1:]
        cov = cov[:r2-1] + str(c2) + cov[r2:]
    cov2 = cov[:]
    cov2 = re.sub("[-]", "", cov2)
    return ("{} [{:>12} : {:>3} restraints touching {} residues]")\
        .format(cov, tbl, len(noe.keys()), len(cov2))


def write_cns_customized_modules(sswt):
    mods = ["scalecoolsetupedited", "scalehotedited"]
    for mod in mods:
        with open(mod) as mod_handle:
            scale_dump = mod_handle.read()
            scale_dump = re.sub("\$sswt", str(sswt), scale_dump)
        os.remove(mod)
        print2file(mod, scale_dump)


def write_cns_gseq(input_file, plen):
    mod = "gseq.inp"
    with open(mod) as gseq_handle:
        gseq_dump = gseq_handle.read()
        gseq_dump = re.sub("\$input_seq2", str(input_file), gseq_dump)
        gseq_dump = re.sub("\$res_renum", str(plen), gseq_dump)
    os.remove(mod)
    print2file(mod, gseq_dump)


def write_cns_dgsa_file(contwt, sswt, mcount, mode,
                        rep1, rep2, mini, f_id, atomselect):
    dgsa_files = ["dgsa.inp"]
    dihed_wt1 = contwt * sswt
    dihed_wt2 = contwt * sswt
    a_select = {1:
                "{===>} md.dg.select=(name ca or name ha or name n or name " +
                "hn\n		        or name c or name cb* or name cg*);",
                2:
                "{===>} md.dg.select=(name ca or name ha or name n or name " +
                "hn\n		        or name c or name cb* or name cg* " +
                "or name o);",
                3:
                "{===>} md.dg.select=(name ca or name ha or name n or name " +
                "hn or name h\n		        or name c or name cb* or " +
                "name cg* or name o);",
                4:
                "{===>} md.dg.select=(name ca or name c or name n or " +
                "name o);",
                5:
                "{===>} md.dg.select=(name ca or name c or name n or name " +
                "o\n		        or name cb);",
                6:
                "{===>} md.dg.select=(name ca or name c or name n or name " +
                "o\n		        or name cb or name h);",
                7:
                "{===>} md.dg.select=(name ca or name c or name n or name " +
                "o\n		        or name cb or name h or name cg*);"
                }
    for dgsa_file in dgsa_files:
        with open(dgsa_file) as dgsa_handle:
            dgsa_dump = dgsa_handle.read()
            dgsa_dump = re.sub("\$contwt", str(contwt), dgsa_dump)
            dgsa_dump = re.sub("\$mcount", str(mcount), dgsa_dump)
            dgsa_dump = re.sub("\$mode", str(mode), dgsa_dump)
            dgsa_dump = re.sub("\$dihed_wt1", str(dihed_wt1), dgsa_dump)
            dgsa_dump = re.sub("\$dihed_wt2", str(dihed_wt2), dgsa_dump)
            dgsa_dump = re.sub("\$rep1", str(int(rep1)), dgsa_dump)
            dgsa_dump = re.sub("\$rep2", str(rep2), dgsa_dump)
            dgsa_dump = re.sub("\$mini", str(mini), dgsa_dump)
            dgsa_dump = re.sub("\$f_id", str(f_id), dgsa_dump)
            dgsa_dump = re.sub("\$a_select", a_select[atomselect], dgsa_dump)
        os.remove(dgsa_file)
        print2file(dgsa_file, dgsa_dump)


def build_models(stage, fasta_file, fasta2_file, ss_file, contwt, sswt, mcount, mode,
                 rep1, rep2, mini, f_id, atomselect, dir_out, cns_suite, debug):
    tbl_list = {}
    for tbl in ["contact.tbl", "ssnoe.tbl", "hbond.tbl", "dihedral.tbl"]:
        if os.path.isfile(tbl):
            tbl_list[tbl] = count_lines(tbl)

    if debug:
        print(seq_fasta(fasta_file))
        if len(fasta2_file)>0:
            print(seq_fasta(fasta2_file))
        if ss_file:
            print(seq_fasta(ss_file))
    flag_dihe = False
    for tbl in sorted(tbl_list.keys()):
        if tbl == "dihedral.tbl":
            flag_dihe = True
        if debug:
            print(coverage_tbl(fasta_file, tbl, flag_dihe,fasta2_file))
        flag_dihe = False

    for filename in glob.glob("iam.*"):
        os.remove(filename)
    write_cns_customized_modules(sswt)
    write_cns_dgsa_file(contwt, sswt, mcount, mode, rep1, rep2,
                        mini, f_id, atomselect)
    for tbl in ["contact.tbl", "ssnoe.tbl", "hbond.tbl", "dihedral.tbl"]:
        if tbl not in tbl_list:
            subprocess.call("sed -i s/{}//g dgsa.inp".format(tbl), shell=True)
    job_file = "#!/bin/bash\n"
    job_file += "echo \"starting cns...\"\n"
    job_file += "touch iam.running\n"
    job_file += "# CNS-CONFIGURATION\n"
    job_file += "source {}/cns_solve_env.sh\n".format(cns_suite)
    job_file += "export KMP_AFFINITY=none\n"
    job_file += "export CNS_CUSTOMMODULE={}/{}\n".format(dir_out, stage)
    job_file += ("{}/intel-x86_64bit-linux/bin/cns_solve < dgsa.inp " +
                 "> dgsa.log \n").format(cns_suite)
    job_file += "if [ -f \"{}_{}.pdb\" ]; then\n".format(f_id, mcount)
    job_file += "   rm iam.running\n"
    job_file += "   echo \"trial structures written.\"\n"
    job_file += "   rm *embed*\n"
    job_file += "   exit\n"
    job_file += "fi\n"
    job_file += "if [ -f \"{}a_{}.pdb\" ]; then \n".format(f_id, mcount)
    job_file += "   rm iam.running\n"
    job_file += "   echo \"accepted structures written.\"\n"
    job_file += "   rm *embed*\n"
    job_file += "   exit\n"
    job_file += "fi\n"
    job_file += "tail -n 30 dgsa.log\n"
    job_file += "echo \"ERROR! Final structures not found!\"\n"
    job_file += "echo \"CNS FAILED!\"\n"
    job_file += "mv iam.running iam.failed\n"
    if os.path.isfile("job.sh"):
        os.remove("job.sh")
    print2file("job.sh", job_file)
    if debug:
        print("Starting job [{}/{}/.job.sh > job.log]".format(dir_out, stage))
    subprocess.call("chmod +x {}/{}/job.sh".format(dir_out, stage),
                    shell=True)
    subprocess.call("{}/{}/job.sh > job.log".format(dir_out, stage),
                    shell=True)
    if os.path.isfile("iam.failed"):
        print("ERROR! Something went wrong while running CNS!")
        print("Check job.log and dgsa.log!")
        sys.exit()


def load_pdb(stage_dir):
    pdb_list = glob.glob(stage_dir + "/*.pdb")
    if not pdb_list:
        pdb_list = glob.glob(stage_dir + "/*.ent")
    if not pdb_list:
        print("ERROR Directory {} has no pdb files!".format(stage_dir))
    return pdb_list


def get_cns_energy(cns_pdb, energy_term):
    if energy_term not in ["overall", "bon", "noe", "vdw"]:
        print("ERROR! Energy term must be one of overall, bon, noe or vdw")
        sys.exit()
    value = "X"
    with open(cns_pdb) as cns_pdb_handle:
        for line in cns_pdb_handle:
            if not line.startswith("REMARK {}".format(energy_term)):
                continue
            value = line.split('=')[1].strip()
            break
    return int(float(value))  # Yes, supposed to only be int, not float


def clash_count(pdb, threshold):
    count = 0
    ca_xyz = xyz_pdb(pdb, "ca")
    for key in sorted([int(x) for x in ca_xyz.keys()]):
        R1 = ca_xyz[str(key)].split()
        x1, y1, z1 = [float(x) for x in R1[:3]]
        for sec_key in sorted([int(x) for x in ca_xyz.keys()]):
            if key >= sec_key:
                continue
            R2 = ca_xyz[str(sec_key)].split()
            x2, y2, z2 = [float(x) for x in R2[:3]]
            d = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
            if d <= threshold:
                count += 1
    return count


def dssp_result(pdb, selection, program_dssp):
    command = '{} {}'.format(program_dssp, pdb) +\
                                ' | grep -C 0 -A 1000 "  #  RESIDUE" ' +\
                                '| tail -n +2'
    process = subprocess.Popen([command],
                               stdout=subprocess.PIPE, shell=True)
    dssp_rows = process.communicate()[0].decode("utf-8")
    RESIDUE = {}
    SS = {}
    PHI = {}
    PSI = {}
    for row in dssp_rows.split('\n'):
        if len(row.strip()) < 1:
            continue
        rnum = row[5:5+6].strip()
        res = row[13].strip()
        sstr = row[16].strip()
        phia = row[103:103+6].strip()
        psia = row[109:109+6].strip()
        res = re.sub("[^A-Z]", "", res)
        sstr = re.sub("[^A-Z]", "", sstr)
        if len(rnum) < 1:
            continue
        rnum = int(rnum)
        if len(res) < 1:
            print("ERROR! Residue not defined for {}".format(res))
            sys.exit()
        if len(phia) < 1:
            print("ERROR! Phi not defined for {}".format(res))
            sys.exit()
        if len(psia) < 1:
            print("ERROR! Psi not defined for {}".format(res))
            sys.exit()
        sstr = "C" if len(sstr) < 1 else sstr
        sstr = re.sub("[\.ISTBG]", "C", sstr)
        RESIDUE[rnum] = res
        SS[rnum] = sstr
        PHI[rnum] = phia
        PSI[rnum] = psia
    if selection == "ss":
        return SS
    elif selection == "phi":
        return PHI
    elif selection == "psi":
        return PSI


def count_ss_match(pdb, fasta_file, ss_file, ss_element, program_dssp):
    if ss_element not in "HEC":
        print("Invalid character! Must be H, E or C")
        sys.exit()
    residue_ss1 = fasta2residues(ss_file)
    residue_ss2 = dssp_result(pdb, "ss", program_dssp)
    count = 0
    for r in residue_ss1.keys():
        if residue_ss1[r] != ss_element:
            continue
        if residue_ss1[r] == residue_ss2[r]:
            count += 1
    return count


def calc_dist(a, b):
    x1, y1, z1 = [float(x) for x in a.split()]
    x2, y2, z2 = [float(x) for x in b.split()]
    return ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5


def ssnoe_tbl_min_pdb_dist(tbl_file, pdb):
    noe_hash = {}
    with open(tbl_file) as tbl_handle:
        for line in tbl_handle:
            row = line.strip()
            row = re.sub("[\)\(]", " ", row)
            C = row.split()
            if C[6] == "or" and C[17] == "or":
                noe_hash[row] = {"left": " ".join([C[2], C[5], C[8], C[11]]),
                                 "right": " ".join([C[13], C[16],
                                                    C[19], C[22]]),
                                 "distance": " ".join([C[23], C[24], C[25]])}
            elif C[6] == "or" and C[17] != "or":
                noe_hash[row] = {"left": " ".join([C[2], C[5], C[8], C[11]]),
                                 "right": " ".join([C[13], C[16]]),
                                 "distance": " ".join([C[17], C[18], C[19]])}
            elif C[6] != "or" and C[11] == "or":
                noe_hash[row] = {"left": " ".join([C[2], C[5]]),
                                 "right": " ".join([C[7], C[10],
                                                    C[13], C[16]]),
                                 "distance": " ".join([C[17], C[18], C[19]])}
            else:
                noe_hash[row] = {"left": " ".join([C[2], C[5]]),
                                 "right": " ".join([C[7], C[10]]),
                                 "distance": " ".join([C[11], C[12], C[13]])}
    xyzPDB = xyz_pdb(pdb, "all")
    for key in noe_hash.keys():
        left = noe_hash[key]["left"]
        right = noe_hash[key]["right"]
        distance = noe_hash[key]["distance"]
        # print(distance)
        L = left.split()
        R = right.split()
        # D = distance.split()
        left_list = {}
        right_list = {}
        for i in range(0, len(L), 2):
            left_list[L[i] + " " + L[i+1]] = 1
        for i in range(0, len(R), 2):
            right_list[R[i] + " " + R[i+1]] = 1
        distance_pdb = 1000.0
        # for k in left_list:
        #     print(k, "=>", left_list[k])
        # sys.exit()
        for le in left_list.keys():
            for ri in right_list.keys():
                L = le.split()
                R = ri.split()
                d = calc_dist(xyzPDB[L[0] + " " +
                              L[1].upper()],
                              xyzPDB[R[0] + " " +
                              R[1].upper()])
                if distance_pdb > d:
                    distance_pdb = d

        noe_hash[key]["pdb_distance"] = distance_pdb
    return noe_hash


def count_satisfied_tbl_rows(pdb, tbl_file, tbl_type, program_dssp):
    if tbl_type not in ["noe", "dihedral"]:
        print("ERROR! Invalid type!")
        sys.exit()
    count = 0
    total = 0
    log_rows = {}
    if tbl_type == "dihedral":
        pass
        # noe = tbl2rows(tbl_file)
        # phi = dssp_result(pdb, "phi", program_dssp)
        # psi = dssp_result(pdb, "psi", program_dssp)
        # for key in noe.keys():
        #     C = key.strip().split()
        #     angle_true = 0.0
        #     if C[5].upper() == "C" and C[10].upper() == "N" and\
        #        C[15].upper() == "CA" and C[20].upper() == "C":
        #         angle_true = float(phi[int(C[2])])
        #     elif C[5].upper() == "N" and C[10].upper() == "CA"\
        #             and C[15].upper() == "C" and C[20].upper() == "N":
        #         angle_true = float(psi[int(C[2])])
        #     else:
        #         print("Undefined dihedral angle")
        #         sys.exit()
        #     viol_flag = 1
        #     d = abs(angle_true - float(C[22]))
        #     d = abs(360.0 - d) if d > 180 else d
        #     if d < (float(C[23]) + 2.0):
        #         count += 1
        #         viol_flag = 0
        #     total += 1

        #     log_rows["{:>3}\t{:.2f}\t{:.2f} # {}".format(viol_flag, d,
        #                                                  angle_true, key)] =\
        #         viol_flag
    else:
        tbl_hash = ssnoe_tbl_min_pdb_dist(tbl_file, pdb)
        for key in sorted(tbl_hash.keys()):
            viol_flag = 1
            distance = tbl_hash[key]["distance"]
            D = [float(x) for x in distance.split()]
            pdb_distance = tbl_hash[key]["pdb_distance"]
            deviation = pdb_distance - (D[0] + D[2])
            if pdb_distance < (D[0] + D[2] + 0.2):
                count += 1
                viol_flag = 0
                deviation = 0.0
            if pdb_distance < (D[0] - D[1] - 0.2):
                count -= 1
                viol_flag = 1
                deviation = -(D[0] - D[1] - pdb_distance)
            log_rows["{:>3}\t{:.2f}\t{:.2f} # {}".format(viol_flag, deviation,
                                                         pdb_distance, key)]\
                = viol_flag
            total += 1
    viol_file = os.path.splitext(os.path.basename(tbl_file))[0] +\
        "_violation.txt"
    write_viol_text = ["#NOE violation check; {} against {}"
                       .format(pdb, tbl_file), "#violation-flag, " +
                       "deviation, actual-measurement, Input-NOE-restraint"]
    for key, _ in sorted(log_rows.items(), key=lambda i: i[0], reverse=True):
        write_viol_text.append(key)
    print2file(viol_file, "\n".join(write_viol_text) + "\n")
    return str(count) + "/" + str(total)


def sum_noe_dev(pdb, tbl_file):
    sum_dev = 0.0
    tbl_hash = ssnoe_tbl_min_pdb_dist(tbl_file, pdb)
    for key in sorted(tbl_hash.keys()):
        viol_flag = 1
        D = [float(x) for x in tbl_hash[key]["distance"].split()]
        pdb_distance = tbl_hash[key]["pdb_distance"]
        if pdb_distance > (D[0] + D[2] + 0.2):
            sum_dev += (pdb_distance - (D[0] + D[2]))
        if pdb_distance < (D[0] - D[1] - 0.2):
            sum_dev += (D[0] - D[1] - pdb_distance)
    return "{:.2f}".format(sum_dev)


def pdb2ss(pdb, program_dssp):
    ss = dssp_result(pdb, "ss", program_dssp)
    ssrow = ""
    for key in sorted(ss.keys()):
        ssrow += ss[key]
    if len(ssrow) < 2:
        print("Looks like DSSP failed!")
        sys.exit()
    return ssrow


def seq_chain(chain):
    seq = ""
    with open(chain) as chain_handle:
        for line in chain_handle:
            if not line.startswith("ATOM"):
                continue
            elif not parse_pdb_row(line.strip(), "aname") == "CA":
                continue
            res = AA3TO1[parse_pdb_row(line.strip(), "rname")]
            seq += res
    return seq


def noe_tbl_violation_coverage(pdb, tbl):
    cov = seq_chain(pdb)
    cov = re.sub("[A-Z]", "-", cov)
    tbl_hash = ssnoe_tbl_min_pdb_dist(tbl, pdb)
    # xyz = xyz_pdb(pdb, "all")
    # print(tbl, pdb)
    for key in sorted(tbl_hash.keys()):
        # print(key)
        left = tbl_hash[key]["left"]
        right = tbl_hash[key]["right"]
        distance = tbl_hash[key]["distance"]
        # print(left)
        # print(right)
        # print(distance)
        L = left.split()
        R = right.split()
        D = distance.split()
        # print(D)
        pdb_distance = tbl_hash[key]["pdb_distance"]
        # print(pdb_distance)
        if pdb_distance > float(D[0]) + float(D[2]) + 0.2:
            cov = cov[:int(L[0]) - 1] + 'x' + cov[int(L[0]):]
        if pdb_distance > float(D[0]) + float(D[2]) + 0.2:
            cov = cov[:int(R[0]) - 1] + 'x' + cov[int(R[0]):]
        if pdb_distance < float(D[0]) - float(D[1]) - 0.2:
            cov = cov[:int(L[0]) - 1] + 'x' + cov[int(L[0]):]
        if pdb_distance < float(D[0]) - float(D[1]) - 0.2:
            cov = cov[:int(R[0]) - 1] + 'x' + cov[int(R[0]):]
    # sys.exit()
    return cov


def assess_dgsa(stage, fasta_file, ss_file, dir_out, mcount, f_id, num_top_models,
                program_dssp, debug):
    seq = seq_fasta(fasta_file)
    pdb_list = load_pdb(os.path.join(dir_out, stage))
    if len(pdb_list) < 2:
        print(("ERROR! Something went wrong while running CNS in {}. Try a " +
               "different atom selection scheme!").format(stage))
        sys.exit()
    if len(pdb_list) < (mcount - 1):
        print(("Warning!! There are some issues! {} models were not " +
               "generated in {}! Try a different atom selection scheme")
              .format(mcount, stage))
    tbl_list = {}
    for tbl in ["contact.tbl", "ssnoe.tbl", "hbond.tbl", "dihedral.tbl"]:
        if os.path.isfile(tbl):
            tbl_list[tbl] = count_lines(tbl)

    for tbl in sorted(tbl_list.keys()):
        if tbl == "dihedral.tbl":
            continue
        search_string = "N1"
        search_string = "N2" if tbl == "ssnoe.tbl" else search_string
        search_string = "HBND" if tbl == "hbond.tbl" else search_string
        if not os.path.isfile("dgsa.log"):
            subprocess.call("touch assess.failed", shell=True)
            print("ERROR! Something went wrong! dgsa.log file not found")
            sys.exit()
        process = subprocess.Popen([('grep NOEPRI dgsa.log | grep {} | ' +
                                     'head -n 1').format(search_string)],
                                   stdout=subprocess.PIPE, shell=True)
        result = process.communicate()[0]
        C = result.split()
        # print(tbl)
        # print(C)
        # print(len(C))
        # print(len(C)-2)
        # print(C[len(C)-2])
        count = int(C[len(C)-2])

        if count != tbl_list[tbl]:
            subprocess.call("touch assess.failed", shell=True)
            print("CNS did not accept all restraints of {}".format(tbl))
            print("Something went wrong, only {}/{} accepted".format(count, count_lines(tbl)))
            sys.exit()
    # Remove "trial" structure of corresponding "accepted" structure
    for i in range(1001):
        if not os.path.isfile("{}a_{}.pdb".format(f_id, i)):
            continue
        else:
            print("Deleting {}_{}.pdb because {}a_{}.pdb exists!"
                  .format(f_id, i, f_id, i))
            os.remove("{}_{}.pdb".format(f_id, i))

    energy_noe = {}
    for pdb in pdb_list:
        if "sub_embed" in pdb or "extended" in pdb:
            continue
        energy_noe[pdb] = get_cns_energy(pdb, "noe")
    # top_pdb = sorted(energy_noe.items(), key=lambda i: [i[1], i[0]])[0][0]

    # if debug:
    #     print("\n\n")
    #     print("           ENERGY            CLASH     SS              NOE " +
    #           "SATISFIED(±0.2A)            SUM OF DEVIATIONS >= 0.2     PDB")
    #     print("--------------------------  -------  -------  ----------------" +
    #           "-----------------------  -------------------------  --------")
    #     print("TOTAL  VDW    BOND   NOE    2.5 3.5  H   E    CONTACTS  SS-NOE" +
    #           "    HBONDS    DIHEDRAL   CONTACTS SS-NOE   HBONDS")

    # for pdb, noe_energy in sorted(energy_noe.items(),
    #                               key=lambda i: i[1]):
    #     e1 = get_cns_energy(pdb, "overall")
    #     e2 = get_cns_energy(pdb, "vdw")
    #     e3 = get_cns_energy(pdb, "bon")
    #     e4 = noe_energy
    #     c1 = clash_count(pdb, 2.5)
    #     c2 = clash_count(pdb, 3.5)
    #     h = count_ss_match(pdb, fasta_file, ss_file, "H", program_dssp)
    #     e = count_ss_match(pdb, fasta_file, ss_file, "E", program_dssp)
    #     if os.path.isfile("contact.tbl"):
    #         n1 = count_satisfied_tbl_rows(pdb, "contact.tbl", "noe", program_dssp)
    #     else:
    #         n1 = "-"
    #     if os.path.isfile("ssnoe.tbl"):
    #         n2 = count_satisfied_tbl_rows(pdb, "ssnoe.tbl", "noe", program_dssp)
    #     else:
    #         n2 = "-"
    #     if os.path.isfile("hbond.tbl"):
    #         n3 = count_satisfied_tbl_rows(pdb, "hbond.tbl", "noe", program_dssp)
    #     else:
    #         n3 = "-"
    #     if os.path.isfile("dihedral.tbl"):
    #         n4 = count_satisfied_tbl_rows(pdb, "dihedral.tbl", "dihedral",
    #                                   program_dssp)
    #     else:
    #         n4 = "-"
    #     if os.path.isfile("contact.tbl"):
    #         s1 = sum_noe_dev(pdb, "contact.tbl")
    #     else:
    #         s1 = "-"
    #     if os.path.isfile("ssnoe.tbl"):
    #         s2 = sum_noe_dev(pdb, "ssnoe.tbl")
    #     else:
    #         s2 = "-"
    #     if os.path.isfile("hbond.tbl"):
    #         s3 = sum_noe_dev(pdb, "hbond.tbl")
    #     else:
    #         s3 = "-"
    #     if debug:
    #         print("{:<6} {:<6} {:<6} {:<6} {:<3} {:<3}  {:<3} {:<3}  "
    #               .format(e1, e2, e3, e4, c1, c2, h, e), end='')
    #         print("{:<9} {:<9} {:<9} {:<9}  {:<8} {:<8} {:<8} {:<25}"
    #               .format(n1, n2, n3, n4, s1, s2, s3, os.path.basename(pdb)))

    # for pdb in sorted(energy_noe.keys(), reverse=True):

    # for tbl in sorted(tbl_list.keys()):
    #     if "dihedral" in tbl:
    #         continue
    #     for pdb in sorted(energy_noe.keys(), reverse=True):
    #         if debug:
    #             print(noe_tbl_violation_coverage(pdb, tbl) + " [ violation of " + os.path.basename(tbl) + " in " + os.path.basename(pdb) + " ]")
    #     if debug:
    #         print()
    i = 1
    return_array = [] 
    for pdb, noe_energy in sorted(energy_noe.items(),
                                  key=lambda i: i[1]):
        if debug:
            print("model{}.pdb <= {}".format(i, pdb))
        shutil.move(pdb, os.path.join(dir_out, "{}".format(stage),"{}_model{}.pdb".format(f_id, i)))
        # shutil.move(pdb, "{}_model{}.pdb".format(f_id, i))
        # print2file(os.path.join(dir_out, "{}_model{}.noe".format(f_id, i)), "{}".format(noe_energy))
        return_array.append(["{}_model{}.noe".format(f_id, i), noe_energy])
        i += 1
        if i > num_top_models:   # Only move the top models
            break
    return return_array


def qa_pcons(dir_out, f_id, program_pcons, debug):
    command = '{} -d {} -A'.format(program_pcons, dir_out)
    process = subprocess.Popen([command],
                               stdout=subprocess.PIPE, shell=True)
    pcons_rows = process.communicate()[0].decode("utf-8")
    return_array = [] 
    for line in pcons_rows.split('\n'):
        if line.startswith(f_id):
            pdb_file, pcons_score = line.split()[:2]
            return_array.append([pdb_file, float(pcons_score)])
    return sorted(return_array, key=lambda x: x[1], reverse=True)


def qa_tmscore(dir_out, native_struct, program_tmscore, debug):
    model_list = load_pdb(dir_out)
    return_array = [] 
    for model in model_list:
        command = '{} {} {}'.format(program_tmscore, model, native_struct)
        process = subprocess.Popen([command],
                                   stdout=subprocess.PIPE, shell=True)
        tmscore_rows = process.communicate()[0].decode("utf-8")
        for line in tmscore_rows.split("\n"):
            if line.startswith("TM-score"):
                tmscore = line.split()[2]
                return_array.append([model.split('/')[-1], float(tmscore)])
    return sorted(return_array, key=lambda x: x[1], reverse=True)
