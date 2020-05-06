import argparse

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", type=str, help="input sequence file in FASTA format")
    parser.add_argument("ss", type=str, help="Secondary structure file in ss2/3 format")
    parser.add_argument("rr", type=str, help="Contact file in CASP format, with or without error as last column")
    parser.add_argument("out_dir", type=str, help="Output directory to write results")
    
    parser.add_argument("-d", "--dist", type=bool, default=False, help="Use predicted distance and error, must use a compatible RR-file")

    parser.add_argument("-om", "--omega", type=str, help="Omega angles file in CASP format")
    parser.add_argument("-th", "--theta", type=str, help="Theta angles file in CASP format")
    parser.add_argument("-noss", "--noss", type=bool, default=False, help="Run without secondary structure")
    parser.add_argument("-ss", "--sstep", type=bool, default=False, help="Save step, on or off")
    parser.add_argument("-tm", "--top_models", type=int, default=5, help="Number of top models selected(from the potentials)")
    parser.add_argument("-rrt", "--rr_type", type=str, default="cb", choices=["ca","cb"], help="Contact type")
    parser.add_argument("-mcount", "--model_count", type=int, default=20, help="Number of potential models generated")
    parser.add_argument("-sr", "--select_rr", type=str, default="all", help="# of contacts from RR-file. Can be 'all' or <float>L")
    parser.add_argument("-lbd", type=float, default=0.4, help="Lambda")
    parser.add_argument("-contwt", type=float, default=10, help="Contact restraint weights")
    parser.add_argument("-sswt", type=float, default=5, help="SS restraint weight")
    parser.add_argument("-rep2", type=float, default=0.85, help="Rep2")
    parser.add_argument("-pthres", type=float, default=7.5, help="Threshold")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = get_args()
