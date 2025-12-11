"""
goal: input: a list of complexes, binary or ternary
(repredict complexes
extract confidence metrics
compute rmsds) in part 2
"""


# import AF3_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from AF3_utils import *
from biopython_utils import *

# load dependancies
import os,  sys
import numpy as np
import argparse
pandas as pd
import time

# general template json file:
template_json_path=f"{SCRIPT_PATH}/no_msa_template.json"
script_start_time = time.time()







######## 

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", nargs="+", type=str, help="List of Input PDBs, used to extract sequences,  then as target template, and as references for RMSDs") # list of pdb files 
parser.add_argument("--target_id", required=True, type=str, help="target id")
parser.add_argument("--smiles", required=False, type=str, help="smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand")
parser.add_argument("--lig_name", required=False, type=str, help="name used for ligand parametrization, eg  FUN ligand")
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file")
#parser.add_argument("--prefix", required=False, type=str, help="prefix for csv file")
args = parser.parse_args()

if ((args.smiles is not None) and (args.lig_name is None)) or ((args.smiles is None) and (args.lig_name is not None)):
  raise ValueError(" if you go for a ternary complex reprediction, please specify smiles AND ligand name")

if (args.smiles is not None) and (args.lig_name is not None):
  smiles=args.smiles
  lig_name=args.lig_name
  ternary=True
  print("reprediction of ternary complexes, ensure the given pdb look like: A target, B binder, ligand")
else:
  ternary=False
  
outdir=args.outdir

# prepare target specific json input file
if ternary:
  target_json=prep_no_msa_target_input_json(gen_template_path=template_json_path, outdir=outdir, target_id=args.target_id, target_cif=args.pdb[0], target_chain_in_cif="A", smiles=smiles, ligand_id=lig_name)
else:
  target_json=prep_no_msa_target_input_json(gen_template_path=template_json_path, outdir=outdir, target_id=args.target_id, target_cif=args.pdb[0], target_chain_in_cif="A", smiles=None, ligand_id=None)
  
##### prepare inputs in the same folder, everything repredicted in one aF3 command
for pdb in args.pdb:
  id=os.path.basename(pdb).replace(".pdb", "")
  structure=load_PDB(pdb)
  seq_binder=chain2seq(structure, "B", make_str=True)
  json=prep_binder_input_json(template_path=target_json, outdir=outdir ,id=id ,seq_binder=seq_binder)


elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
n_binder=len(args.pdb)
print(f"Finished AF3 input preparation for {n_binder} complexes. Script execution took: "+elapsed_text)


