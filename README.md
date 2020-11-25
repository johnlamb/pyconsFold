# pyconsFold:
A python modelling framework built on top of CNS.
Support for both trRosetta distance predictions and CASP format contact predictions, both binary and distance based.

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
