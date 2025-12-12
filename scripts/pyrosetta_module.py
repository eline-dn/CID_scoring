# load dependencies
import os,  sys
import numpy as np
import argparse
import pandas as pd
import time
# import pyrosetta_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from pyrosetta_utils import *

script_start_time = time.time()
# one binder scoring function:
def score_target_n_binder(pdb):
  data={}
  pose=load_pose(pdb, assert_lig=False)
  # interchain contacts computations:
  #data["atom_contact_BT1"]=count_chain_atom_contacts(pose=pose, ch_id1=1, ch_id2=2, delta=0.2)
  data["atom_contact_BT2"]= get_atomic_contact_data(pose)
  data["atom_contact_BT3"]= count_atomic_contacts(pdb=pdb, chain1="A", chain2="B", delta=0.2)
  interface_scores, interface_AA, interface_residues_pdb_ids_str=score_nolig_interface(pdb_file=pdb, interface= "A_B", binder_chain="B")
  data={**data, **interface_scores}
  return data

def score_ternary_complex(pdb):
  data={}
  pose=load_pose(pdb, assert_lig=True) # load pose and check that the last residue in the pose is a ligand
  # derive binary pose from the ternary complex: (i.e. ligand + binder)
  binder_start = pose.conformation().chain_begin(2)
  lb_pose = pyrosetta.rosetta.core.pose.Pose()
  pyrosetta.rosetta.core.pose.append_subpose_to_pose(lb_pose, pose, binder_start, pose.size(), 1)
  
  # interchain contacts computations:
  #data["atom_contact_BL1"]=count_chain_atom_contacts(pose=pose, ch_id1=2, ch_id2=3, delta=0.2)
  data["atom_contact_BL2"]= get_atomic_contact_data(lb_pose) # this one is probably going to separate both and give zero, do not lend it to much confidence
  data["atom_contact_BL3"]= count_atomic_contacts(pdb=pdb, chain1="B", chain2="L", delta=0.2)
  interface_scores, interface_AA, interface_residues_pdb_ids_str=score_lig_interface(pdb_file, interface= "AF_B", binder_chain="B")
  data={**data, **interface_scores}
  return data  


######## 

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", nargs="+", type=str, help="List of Input PDBs") # list of pdb files 
parser.add_argument("--mk_params", type=bool, help=" need to create a Params files in the pdb files'dir? Set to true if the pdbs are ternary complexes. in that case, add a smiles for the ligand")
parser.add_argument("--smiles", required=False, type=str, help="smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand")
parser.add_argument("--lig_name", required=False, type=str, help="name used for ligand parametrization, eg  FUN ligand")
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file")
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

#  init pyrosetta and parametrize ligand (or not):
if args.mk_params:
  params=create_param(smiles=args.smiles, pdb_file=args.pdb[0], lig_name=args.lig_name)
  pyr_init(params=params)
else: 
  pyr_init() # initialise without the param file

for pdb in args.pdb:
  id=os.path.basename(pdb).replace(".pdb", "")
  if args.mk_params:
    data=score_ternary_complex(pdb)
  else: 
    data=score_target_n_binder(pdb)
  # write to csv
  data["id"]=id
  df=pd.DataFrame(data=data, index=[data["id"]])
  csv_path=f"{outdir}/{prefix}_pyrosetta_metrics.csv"
  df.to_csv(csv_path, mode="a", index=False, header=not pd.io.common.file_exists(csv_path))
    
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
n_binder=len(args.pdb)
print(f"Finished Pyrosetta scoring for {n_binder} complexes. Script execution took: "+elapsed_text)
