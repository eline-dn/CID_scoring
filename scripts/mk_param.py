import os, sys, glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.PDB
import json


# import AF3_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from AF3_utils import *
from biopython_utils import *

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--lig_name', type=str, required=True, help='ligand name')
parser.add_argument('--smiles', type=str, required=True, help='ligand smiles as a string')
parser.add_argument('--cif_path', type=str, required=True, help='template cif file (output from af3) to use for atom numbering')
args = parser.parse_args()



pdb_file=af3_out_2_norm_pdb(args.cif_path, lig=True, lig_name=None)


# param file:

# trying to generate a param file from a pdb with ligand and smiles. goal: atoms names match the pdb names, (else pyrosetta crashes)

import pyrosetta
pyrosetta.init(extra_options='-mute all') # required for test
from rdkit_to_params import Params
# eg: "c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl"
#smiles="c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl"
p = Params.from_smiles_w_pdbfile(pdb_file, args.smiles, name=args.lig_name) # the name has to match.
print(p.is_aminoacid()) # False
p.dump(f'{args.lig_name}.params')
p.test().dump_pdb('test.pdb')
print("Check the ligand atom numbering in test.pdb, the param file and the af3 outputs: they have to match")
