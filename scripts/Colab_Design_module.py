"""
goal: from reference pdbs, repredict the complex without the ligand and retrive the reprediction metrics
"""
# load dependancies
import os, sys
import numpy as np
import argparse
import pandas as pd
import time
import glob

# import AF3_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from Colab_utils import *
from biopython_utils import *


######## 
script_start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("--ref_pdb", nargs="+", type=str, help="List of ref PDBs, used to extract sequences, target template, and as references for RMSDs. Ensure the chain id's are correct.") 
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file and repredictions")
parser.add_argument("--prefix", required=False, type=str, help="prefix for csv file")
args = parser.parse_args()

  
if args.outdir is not None:
  outdir=args.outdir
else: 
  outdir="."

if args.prefix is not None:
  prefix=args.prefix
else: 
  prefix=""


for pdb in args.ref_pdb:
  id=os.path.basename(pdb).replace(".pdb", "")
  # extract template for the target:
  structure=load_PDB(pdb)
  template_path=f"{outdir}/{id}_target_template.pdb"
  extract_template_target(structure, template_path)
  # binder sequence and length:
  binder_length=chain2length(structure, "B")
  binder_sequence=chain2seq(structure,  "B", True)
  # run reprediction:
  model=init_colab()
  data=run_prediction(model, template_path, binder_length, binder_sequence, outdir, id)
  
  # compute RMSDs
  rmsds={}
  for i, model in enumerate(glob.glob(f"{outdir}/{id}_model*.pdb")):
    ref= load_PDB(pdb)
    mov=load_PDB(model)
    # chain structure
    mapping={"A":"A", "B":"B"}
    # compute RMSDs
    rmsds[i]=binary_RMSDs(ref, mov, mapping)
    
  # mean of the two models:
  rmsd={}
  for key in rmsds[0].keys():
    rmsd[key]=round((rmsds[0][key] + rmsds[1][key])/2, 2)
  print(data)
  print(rmsd)
  # save in csv
  df=pd.DataFrame(data={**data, **rmsd}, index=[data["id"]])
  csv_path=f"{outdir}/{prefix}_Colab_Design_reprediction_metrics.csv"
  df.to_csv(csv_path, mode="a", index=False, header=not pd.io.common.file_exists(csv_path))
    
  


elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
n_binder=len(args.ref_pdb)
print(f"Finished  Colab Design structure reprediction output scoring for {n_binder} complexes. Script execution took: "+elapsed_text)
