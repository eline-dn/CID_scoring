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
import glob

# import AF3_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from AF3_utils import *
from biopython_utils import *


######## 
script_start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("--ref_dir", required=True, type=str, help="Directory containing reference PDB or CIF files.") 
parser.add_argument("--binary_out_dir", required=True, type=str, help="Directory containing binary AF3 output folders.")
parser.add_argument("--ternary_out_dir", required=True, type=str, help="Directory containing ternary AF3 output folders.")
parser.add_argument("--data_csv", required=True, type=str, help="Path to data_complexes.csv file with chain information.")
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file")
parser.add_argument("--prefix", required=False, type=str, help="prefix for csv file")
args = parser.parse_args()

# Usage example:
# python AF3_output.py --binary_out_dir /work/lpdi/users/eline/CID_scoring/benchmark/binary_out --ternary_out_dir /work/lpdi/users/eline/CID_scoring/benchmark/ternary_out --ref_dir /work/lpdi/users/eline/CID_scoring/benchmark/pdb_files --data_csv /work/lpdi/users/eline/CID_scoring/benchmark/data_complexes.csv --outdir /work/lpdi/users/eline/CID_scoring/benchmark --prefix ""

# Load the data CSV
data_df = pd.read_csv(args.data_csv)
data_df['pdb_id'] = data_df['pdb_id'].str.upper()  # Ensure uppercase for matching
  
if args.outdir is not None:
  outdir=args.outdir
else: 
  outdir="."

if args.prefix is not None:
  prefix=args.prefix
else: 
  prefix=""

# Get all IDs from binary and ternary directories
binary_ids = set()
if os.path.exists(args.binary_out_dir):
    binary_ids = {d[:4].upper() for d in os.listdir(args.binary_out_dir) if os.path.isdir(os.path.join(args.binary_out_dir, d))}

ternary_ids = set()
if os.path.exists(args.ternary_out_dir):
    ternary_ids = {d[:4].upper() for d in os.listdir(args.ternary_out_dir) if os.path.isdir(os.path.join(args.ternary_out_dir, d))}

all_ids = binary_ids.union(ternary_ids)

all_data = []

for id in all_ids:
  print("Processing id:", id)
  # Get chain info from CSV
  row = data_df[data_df['pdb_id'].str.upper() == id.upper()]
  if row.empty:
    print(f"No data found for id {id} in CSV, skipping")
    continue
  chain_a = row['chain_a'].values[0]
  chain_b = row['chain_b'].values[0]
  ligand_id = row.iloc[0, 1]  # Assuming ligand_id is the second column
  # find reference:
  ref_file = None
  for case_id in [id.upper(), id.lower()]:
    for ext in ['.pdb', '.cif']:
      potential_file = os.path.join(args.ref_dir, f"{case_id}{ext}")
      if os.path.exists(potential_file):
        ref_file = potential_file
        break
    if ref_file:
      break
  if ref_file is None:
    print(f"no reference file for id {id}, continuing")
    continue
  
  # load the reference structure:
  try:
    if ref_file.endswith('.pdb'):
      ref = load_PDB(ref_file)
    elif ref_file.endswith('.cif'):
      ref = load_CIF(ref_file)
    else:
      raise ValueError(f"Unsupported reference file format for {ref_file}")
  except Exception as e:
    print(f"Error loading reference for {id}: {e}, skipping")
    continue
  
  # Process binary if exists
  if id in binary_ids:
    try:
      # Find the full dir name
      full_dirs = [d for d in os.listdir(args.binary_out_dir) if d.upper().startswith(id.upper()) and os.path.isdir(os.path.join(args.binary_out_dir, d))]
      if not full_dirs:
        raise ValueError(f"No binary dir found for {id}")
      dir_path = os.path.join(args.binary_out_dir, full_dirs[0])
      dir_name = os.path.basename(dir_path)
      confidence = os.path.join(dir_path, f"{dir_name}_summary_confidences.json")
      cif = os.path.join(dir_path, f"{dir_name}_model.cif")
      if not os.path.exists(confidence) or not os.path.exists(cif):
        raise ValueError(f"Missing files for binary {id}")
      data = extract_af3_confidence_metrics(confidence)
      movpdb = af3_out_2_norm_pdb(cif, lig=False, lig_name=None)
      mov = load_PDB(movpdb)
      mapping = {chain_a: "A", chain_b: "B"}
      rmsds = binary_RMSDs(ref, mov, mapping)
      # add prefix
      data = {f"bin_{k}": v for k, v in data.items()}
      rmsds = {f"bin_{k}": v for k, v in rmsds.items()}
      # contacts
      binder_contacts = hotspot_residues(movpdb, binder_chain="B", atom_distance_cutoff=4.0)
      binder_contacts_n = len(binder_contacts)
      data["bin_n_contacts_target_binder"] = binder_contacts_n
      data["id"] = id
      data["type"] = "binary"
      data["error"] = ""
      data = {**data, **rmsds}
      all_data.append(data)
    except Exception as e:
      print(f"Error processing binary for {id}: {e}")
      # Create failed entry
      data = {"id": id, "type": "binary", "error": str(e)}
      # Add NaN for all metrics
      # Assuming standard keys, but to be general, we can list them
      metric_keys = ["bin_iptm", "bin_ptm", "bin_mean_plddt", "bin_target_plddt", "bin_binder_plddt", "bin_target_iptm", "bin_binder_iptm", "bin_AA_pair_iptm", "bin_AB_pair_iptm", "bin_BB_pair_iptm", "bin_AB_pair_pae_min", "bin_BA_pair_pae_min", "bin_target_RMSD", "bin_binder_RMSD", "bin_full_RMSD", "bin_d_binding_site", "bin_n_contacts_target_binder"]
      for key in metric_keys:
        data[key] = np.nan
      all_data.append(data)
  
  # Process ternary if exists
  if id in ternary_ids:
    try:
      # Find the full dir name
      full_dirs = [d for d in os.listdir(args.ternary_out_dir) if d.upper().startswith(id.upper()) and os.path.isdir(os.path.join(args.ternary_out_dir, d))]
      if not full_dirs:
        raise ValueError(f"No ternary dir found for {id}")
      dir_path = os.path.join(args.ternary_out_dir, full_dirs[0])
      dir_name = os.path.basename(dir_path)
      confidence = os.path.join(dir_path, f"{dir_name}_summary_confidences.json")
      cif = os.path.join(dir_path, f"{dir_name}_model.cif")
      if not os.path.exists(confidence) or not os.path.exists(cif):
        raise ValueError(f"Missing files for ternary {id}")
      # Detect ligand chain
      struct = load_CIF(cif)
      chains = [chain.id for chain in struct[0]]
      lig_chains = [c for c in chains if c not in ['A', 'B']]
      if len(lig_chains) != 1:
        raise ValueError(f"Could not detect single ligand chain for {id}, found {lig_chains}")
      lig_name = lig_chains[0]
      # Find ligand chain in reference
      ref_lig_chain = None
      for chain in ref[0]:
        for res in chain:
          if res.resname.strip() == ligand_id:
            ref_lig_chain = chain.id
            break
        if ref_lig_chain:
          break
      if not ref_lig_chain:
        raise ValueError(f"Could not find ligand chain in reference for {id} with ligand_id {ligand_id}")
      data = extract_af3_confidence_metrics(confidence, lig_name=lig_name)
      movpdb = af3_out_2_norm_pdb(cif, lig=True, lig_name=lig_name)
      mov = load_PDB(movpdb)
      mapping = {chain_a: "A", chain_b: "B", ref_lig_chain: "L"}
      rmsds = ternary_RMSDs(ref, mov, mapping)
      if rmsds["ligand_rmsd"] == "failed":
        raise ValueError(f"Bad ligand reprediction for {id}")
      # add prefix
      data = {f"ter_{k}": v for k, v in data.items()}
      rmsds = {f"ter_{k}": v for k, v in rmsds.items()}
      # contacts
      binder_contacts = hotspot_residues(movpdb, binder_chain="B", atom_distance_cutoff=4.0)
      binder_contacts_n = len(binder_contacts)
      data["ter_n_contacts_target_binder"] = binder_contacts_n
      data["id"] = id
      data["type"] = "ternary"
      data["error"] = ""
      data = {**data, **rmsds}
      all_data.append(data)
    except Exception as e:
      print(f"Error processing ternary for {id}: {e}")
      # Create failed entry
      data = {"id": id, "type": "ternary", "error": str(e)}
      # Add NaN for all metrics
      metric_keys = ["ter_iptm", "ter_ptm", "ter_mean_plddt", "ter_target_plddt", "ter_binder_plddt", "ter_ligand_plddt", "ter_target_iptm", "ter_binder_iptm", "ter_ligand_iptm", "ter_AA_pair_iptm", "ter_AB_pair_iptm", "ter_BB_pair_iptm", "ter_AL_pair_iptm", "ter_BL_pair_iptm", "ter_LL_pair_iptm", "ter_AB_pair_pae_min", "ter_BA_pair_pae_min", "ter_AL_pair_pae_min", "ter_LA_pair_pae_min", "ter_BL_pair_pae_min", "ter_LB_pair_pae_min", "ter_target_RMSD", "ter_binder_RMSD", "ter_full_RMSD", "ter_d_binding_site", "ter_ligand_rmsd", "ter_n_contacts_target_binder"]
      for key in metric_keys:
        data[key] = np.nan
      all_data.append(data)

# Save to CSV
if all_data:
  df = pd.DataFrame(all_data)
  csv_path = f"{outdir}/{prefix}_AF3_metrics.csv"
  df.to_csv(csv_path, index=False)
  print(f"Saved metrics to {csv_path}")
else:
  print("No data to save")

