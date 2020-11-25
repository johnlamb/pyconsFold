# pyconsFold:
A python modelling framework built on top of CNS.
Support for both trRosetta distance predictions and CASP format contact predictions, both binary and distance based.

## Install
pyconsFold require a working installation of [CNS](http://cns-online.org). This needs to be done manually due to license.

1. Install CNS
	1. [Request](http://cns-online.org/cns_request/) a download link from CNS.
	2. Follow the emailed instructions to download "cns_solve_1.3_all_intel-mac_linux.tar.gz
	3. Extract the files `tar xzvf cns_solve_1.3_all_intel-mac_linux.tar.gz`
	4. Change into the resulting directory `cd cns_solve_1.3`
	5. Unhide the bash-specific file `mv .cns_solve_env_sh cns_solve_env.sh`
	6. In this resulting file, replace `_CNSsolve_location_` with the CNS installation folder. If you extracted the file in your homefolder then the CNS installation would be: `/home/<your username>/cns_solve_1.3`
	7. Source CNS, `source cns_solve_env.sh`, to make this permanent and to prevent you having to do this every time, add it to your .bashrc file.
	8. Test CNS by going into the test folder `cd test` and run the tests `../bin/run_tests -tidy *.inp`

2. Install pyconsFold
	1. Run `pip install pyconsFold`

3. Optional: If you clone this github repo, you can run a suite of tests using `python3 run_test.py`

## Usage
```python
import pyconsFold

pyconsFold.model(fasta, contacts, out_dir)
```

## Functions
	model      -- Classic modelling using binary contact predictions (although the contact file can contain distance and errors they wont be used)
	model_dist -- Model using distance and errors, requires either a CASP-formated rr file with additional column with standard error in Ångströms or a trRosetta-contacts file in npz-format.
	model_dock -- Perform modelling and docking of two protein chains. Requires _one_ contacts file with both inter- and intra-contacts.

## Utilities
QA-function arguments to all above functions:
	pcons (default False) -- If set to true, gives pcons score for all models (using either pcons installed in the PATH or the builtin binary)
	tmscore_pdb_file      -- If a structure file is supplied, runs all models against this (presumed) native structure and reports the TMscore (using either TMscore in the PATH or builtin binary)

## Extras
```python
from pyconsFold.utils import npz_to_casp, pdb_to_npz

npz_to_casp("trRosetta.npz")  ##  Converts trRosetta distance and angle predictions to CASP format in separate files

pdb_to_npz("structure.pdb")   ##  Converts a structure (pdb/mmCif) to trRosetta distances and angles, useful when investigating how well a model conforms to restraints
```

To install: `pip install -e .`
