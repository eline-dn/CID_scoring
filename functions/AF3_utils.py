import os, sys, glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import getpass
import subprocess
import time
import importlib
from shutil import copy2
import json
import ast
# add utils to path
from biopython_utils import *

"""
prep json files from seq, template, smiles
=> run from a dir of json files
extract relevant metrics from json output file
-save what is needed for ipSAE/ LIS 
"""

def prep_binder_input_json(template_path, outdir,id,seq_binder):
  with open(template_path, "r") as f:
    data = json.load(f)
  output_path = f"{outdir}/{id}.json" # binder sequence specific json input for each af3 reprediction
  data["name"]=id
  for entry in data["sequences"]:
      if "protein" in entry:
          if entry["protein"]["id"] == ["B"]:
              entry["protein"]["sequence"] = seq_binder

  # write modified JSON to a new file
  with open(output_path, "w") as f:
      json.dump(data, f, indent=2)

  print(f"json input file saved in {output_path}")
  return(output_path)

def prep_no_msa_target_input_json(gen_template_path, outdir, target_id, target_cif, target_chain_in_cif="A", smiles=None, ligand_id=None):
  with open(gen_template_path, "r") as f:
    data = json.load(f)
  output_path = f"{outdir}/{target_id}.json" # target sequence specific json input, to be modified for each binder
  data["name"]=target_id
  for entry in data["sequences"]:
      if "protein" in entry:
          if entry["protein"]["id"] == ["A"]:
            struct=load_CIF(target_cif)
            target_seq=chain2seq(struct, chain_id=target_chain_in_cif, make_str=True)
            entry["protein"]["sequence"] = target_seq
            target_len=chain2length(struct,  chain_id=target_chain_in_cif)
            # we assume every template has target in chain A, correctly numbered and same length as known seq
            Indices=list(range(target_len))
            target_template_dict={"mmcifPath":target_cif,
                                    "queryIndices":Indices ,
                                    "templateIndices":Indices}
              
            entry["protein"]["templates"] = [target_template_dict]
  if (ligand_id is not None) and (smiles is not None):
    dict={"ligand":{"id":ligand_id,
                                 "smiles":smiles}}
    data["sequences"].append(dict)
  
  # write modified JSON to a new file
  with open(output_path, "w") as f:
      json.dump(data, f, indent=2)
  print(f"json input file template saved in {output_path}")
  return(output_path)

def extract_af3_confidence_metrics(confidence, lig_name=None):
  # from confidence *_summary_confidences.json file
  #load JSON
  with open(confidence, "r") as f:
      data = json.load(f) # data["chain_iptm", "chain_pair_iptm", "chain_pair_pae_min","chain_ptm","fraction_disordered": 0.0, "has_clash": 0.0, "iptm": 0.7, "ptm": 0.75, "ranking_score": 0.71]
  id=os.path.basename(confidence)
  id=id.replace("_summary_confidences", "") # e.g. t2_2_7_4_t0.3_ltp0.2_dcut500.0
  id=id.replace(".json", "")
  data["id"]=id
  # also get atoms plddts from ..._confidences.json:
  conf=confidence.replace("_summary", "")
  with open(conf,"r") as f:
      conf=json.load(f)
  #data["atom_plddt"]=conf["atom_plddts"]
  data["mean_plddt"]=sum(conf["atom_plddts"]) / len(conf["atom_plddts"])+0.00001
  #data["atom_chain_ids"]=conf["atom_chain_ids"]
  data["chain_plddt"]={}
  chain_atom_len={}
  for i,chain_id in enumerate(conf["atom_chain_ids"]):
      if chain_id not in data["chain_plddt"].keys():
          data["chain_plddt"][chain_id]=0
      data["chain_plddt"][chain_id]+=conf["atom_plddts"][i]
      if chain_id not in chain_atom_len.keys():
          chain_atom_len[chain_id]=0
      chain_atom_len[chain_id]+=1
  for chain, plddt in data["chain_plddt"].items():
      data["chain_plddt"][chain]=plddt/chain_atom_len[chain] # = mean atom plddt per chain
  ## split col chain_plddt into binder_plddt, ligand_plddt; target_plddt
  data["chain_plddt"] = data["chain_plddt"].apply(
    lambda x: ast.literal_eval(x) if isinstance(x, str) else x
  )
  data["binder_plddt"] = data["chain_plddt"].apply(lambda d: d.get("B"))
  if lig_name is not None:
    data["ligand_plddt"] = data["chain_plddt"].apply(lambda d: d.get(lig_name))
  data["target_plddt"] = data["chain_plddt"].apply(lambda d: d.get("A"))
  del data["chain_plddt"]
  # same for chain_iptm and chain_ptm
  #for 
  return(data)

