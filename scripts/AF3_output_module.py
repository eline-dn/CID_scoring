"""
(goal previous scripts: input: a list of complexes, binary or ternary
repredict complexes)
what is done here:
extract confidence metrics
compute rmsds
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

###### a few helper functions


######## 
script_start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("--ref_pdb", nargs="+", type=str, help="List of ref PDBs, used to extract sequences, as target template, and as references for RMSDs. Ensure the id is the same as the folder with repredicted strucures and confidences. Also ensure the chain id's are matching.") 
parser.add_argument("--AF3_outs", nargs="+", type=str, help="List of AF3 output folders. Ensure the id is the same as the ref_pdb.") 
#parser.add_argument("--target_id", required=True, type=str, help="target id")
#parser.add_argument("--smiles", required=False, type=str, help="smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand")
parser.add_argument("--lig_name", required=False, type=str, help="name used for ligand parametrization, eg  FUN ligand")
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file")
#parser.add_argument("--prefix", required=False, type=str, help="prefix for csv file")
args = parser.parse_args()

if (args.lig_name is not None):
  lig_name=args.lig_name
  ternary=True
  print("reprediction of ternary complexes, ensure the given pdb look like: A target, B binder, ligand")
  metrics_prefix="ter_"
else:
  ternary=False
  metrics_prefix="bin_"
  
outdir=args.outdir


##### prepare inputs in the same folder, everything repredicted in one aF3 command
for pdb in args.pdb:
  id=os.path.basename(pdb).replace(".pdb", "")
  structure=load_PDB(pdb)
  seq_binder=chain2seq(structure, "B", make_str=True)
  json=prep_binder_input_json(template_path=target_json, outdir=outdir ,id=id ,seq_binder=seq_binder)


elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
n_binder=len(args.pdb)
print(f"Finished structure reprediction output scoring for {n_binder} complexes. Script execution took: "+elapsed_text)

