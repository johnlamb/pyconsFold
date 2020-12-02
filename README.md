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
	7. Source CNS, `source cns_solve_env.sh`\*, to make this permanent and to prevent you having to do this every time, add it to your .bashrc file.
	8. Test CNS by going into the test folder `cd test` and run the tests `../bin/run_tests -tidy *.inp`
_* If you get an error about csh interpreter, you need to install csh_

2. Install pyconsFold
	1. Run `pip install pyconsFold`

3. Optional: If you clone this github repo, you can run a suite of tests using `python3 run_test.py`

## Usage
```python
import pyconsFold

pyconsFold.model_dist(fasta, contacts, out_dir)
```

## Functions
	model      -- Classic modelling using binary contact predictions (although the contact file can contain distance and errors they wont be used)
	model_dist -- Model using distance and errors, requires either a CASP-formated rr file with additional column with standard error in Ångströms or a trRosetta-contacts file in npz-format.
	model_dock -- Perform modelling and docking of two protein chains. Requires _one_ contacts file with both inter- and intra-contacts.
### Top arguments
```
rr_pthres	--	Threshold for the confidence we want in a prediction (default model(0.80), model_dist(0.45), model_dock(0.50))
rr_sep		--	Separation between contacts (default 0)
save_step	--	Save working steps (default False)
stage2		--	Run stage2, filter contacts vs generated structure and generate new structures with filtered contacts (default False)
debug		--	Write out debug information (default False)
selectrr	--	How many contacts to use? Can be "all", "#L", or #. (default "all")
mcount		--	How many models to generate? (default 20)
top_models	--	How many of the generated models should be ranked and saved? (default 20)
use_angles	--	If predicted angels should be used, only works with npz (default False)
omega		--	RR-formated file with omega angles (if npz are not used) (default '')
theta		--	RR-formated file with theta angles (if npz are not used) (default '')
```

## Utilities
QA-function arguments to all above functions:
* pcons (default False) -- If set to true, gives pcons score for all models (using either pcons installed in the PATH or the builtin binary)
* tmscore_pdb_file      -- If a structure file is supplied, runs all models against this (presumed) native structure and reports the TMscore (using either TMscore in the PATH or builtin binary)

## Extras
```python
from pyconsFold.utils import npz_to_casp, pdb_to_npz

npz_to_casp("trRosetta.npz")  ##  Converts trRosetta distance and angle predictions to CASP format in separate files

pdb_to_npz("structure.pdb")   ##  Converts a structure (pdb/mmCif) to trRosetta distances and angles, useful when investigating how well a model conforms to restraints
```

## Adjustable parameters for CNS, advanced
```
rrtype		--	Between which atoms in a residue are the contacts? (default 'cb')
lbd			--	Lambda, 0.1-10 (default 0.4)
contwt		--	Contact restraint weights, 0.1-10000 (default 10)
sswt		--	Secondary structure weights, 0.1-100 (default 5)
bin_values	--	Dictionary of bin_values for converstion of npz to RR-format, see source code (default {})



```
