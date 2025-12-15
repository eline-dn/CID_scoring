"""
(goal previous scripts: input: a list of complexes, binary or ternary
repredict complexes)
what is done here:
extract confidence metrics
compute rmsds
"""
# load dependencies
import os,  sys
import numpy as np
import argparse
import pandas as pd
import time

# import AF3_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from AF3_utils import *
from biopython_utils import *


######## 
script_start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("--ref_pdb", nargs="+", type=str, help="List of ref PDBs, used to extract sequences, as target template, and as references for RMSDs. Ensure the id is the same as the folder with repredicted strucures and confidences. Also ensure the chain id's are matching.") 
parser.add_argument("--AF3_outs", nargs="+", type=str, help="List of AF3 output folders. Ensure the id is the same as the ref_pdb.") 
#parser.add_argument("--target_id", required=True, type=str, help="target id")
#parser.add_argument("--smiles", required=False, type=str, help="smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand")
parser.add_argument("--lig_name", required=False, type=str, help="name used for ligand parametrization, eg  FUN ligand")
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file")
parser.add_argument("--prefix", required=False, type=str, help="prefix for csv file")
args = parser.parse_args()

if (args.lig_name is not None):
  lig_name=args.lig_name
  ternary=True
  print("reprediction of ternary complexes, ensure the given pdb look like: A target, B binder, ligand")
  metrics_prefix="ter_"
else:
  ternary=False
  metrics_prefix="bin_"
  
if args.outdir is not None:
  outdir=args.outdir
else: 
  outdir="."

if args.prefix is not None:
  prefix=args.prefix
else: 
  prefix=""


for dir in args.AF3_outs:
  id=dir.split("/")[-1]
  # extract confidence metrics:
  confidence=dir + "/" +id+"_summary_confidences.json"
  data=extract_af3_confidence_metrics(confidence)
  # find reference:
  for pdb in args.ref_pdb:
    if id in pdb:
      ref_pdb_file=pdb
  # load the structures:
  ref= load_PDB(ref_pdb_file)
  # TO DO: sanitize structure fiel with
  cif=os.path.join(dir,f"{id}_model.cif")
  # chain structure
  mapping={"A":"A", "B":"B", "L":args.lig_name}
  # compute RMSDs
  if ternary:
    movpdb=af3_out_2_norm_pdb(cif, lig=True, lig_name=args.lig_name)
    mov= load_PDB(movpdb)
    rmsds=ternary_RMSDs(ref, mov, mapping)
  else:
    movpdb=af3_out_2_norm_pdb(cif, lig=False, lig_name=None)
    mov= load_PDB(movpdb)
    rmsds=binary_RMSDs(ref, mov, mapping)
  # add a prefix on the metrics names:
  data = {f"{metrics_prefix}{k}": v for k, v in data.items()}
  rmsds = {f"{metrics_prefix}{k}": v for k, v in rmsds.items()}
  # add the id:
  data["id"]=id
  # save in csv
  df=pd.df([{**data, **rmsds}])
  csv_path=f"{outdir}/{prefix}_AF3_{metrics_prefix}reprediction_metrics.csv"
  df.to_csv(csv_path, mode="a", index=False, header=not pd.io.common.file_exists(csv_path))
    
  


elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
n_binder=len(args.AF3_outs)
print(f"Finished structure reprediction output scoring for {n_binder} complexes. Script execution took: "+elapsed_text)

