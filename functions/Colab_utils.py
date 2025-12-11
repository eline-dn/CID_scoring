


# import dependencies:
import pandas as pd
from Bio.PDB.Polypeptide import is_aa
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
import Bio.PDB
from scipy.special import softmax
from colabdesign import mk_afdesign_model, clear_mem
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.loss import get_ptm, mask_loss, get_dgram_bins, _get_con_loss, get_plddt_loss, get_exp_res_loss, get_pae_loss, get_con_loss, get_rmsd_loss, get_dgram_loss, get_fape_loss
from colabdesign.shared.utils import copy_dict
from colabdesign.shared.prep import prep_pos
from Bio.PDB import PDBIO, StructureBuilder, PDBParser

# source functions:
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from biopython_utils import *


def extract_template_target(structure, output_pdb_path):
  target= copy_structure_with_only_chain(structure, "A")
  write_pdb(structure, output_pdb_path)
  print(f"target template saved in {output_pdb_path}")


def init_colab():
  clear_mem()
  params = '/work/lpdi/users/goldbach/software/colabdesign/params' 
  #run reprediction with Colab Design:
  # compile complex prediction model
  model = mk_afdesign_model(protocol="binder", num_recycles=3, data_dir=params, 
                                              use_multimer=True,
                                              use_templates=True,
                                               use_initial_guess=False, #Introduce bias by providing binder atom positions as a starting point for prediction.
                                               use_initial_atom_pos=False) # Introduce atom position bias into the structure module for atom initilisation.
  return(model)

def run_prediction(model, template, binder_length, binder_sequence, outdir, id):
  model.prep_inputs(pdb_filename=template,
                        chain="A",
                        #binder_chain="B",# do not specifiy if the template only contains the target
                        binder_len=binder_length,
                        rm_target_seq=False, #b
                        use_binder_template=False #a
                        #,rm_template_ic=True #c
                        )
  
  prediction_stats = {}
  for model_num in [0,1]:
    model.predict(seq=binder_sequence,
                  models=[model_num],
                  num_recycles=3)

    predicted_complex_pdb = f"{outdir}/{id}_model{model_num}.pdb"
    model.save_pdb(predicted_complex_pdb)
    prediction_metrics = copy_dict(model.aux["log"]) # contains plddt, ptm, i_ptm, pae, i_pae

    # extract the statistics for the model
    stats = {
        f"pLDDT": round(prediction_metrics['plddt'], 2),
        f"pTM": round(prediction_metrics['ptm'], 2),
        f"i_pTM": round(prediction_metrics['i_ptm'], 2),
        f"pAE": round(prediction_metrics['pae'], 2),
        f"i_pAE": round(prediction_metrics['i_pae'], 2)
          }       


    prediction_stats[model_num+1] = stats # 2 dictionnaries index 1 and 2 to eventually add to the metrics df

  data={}
  for key in prediction_stats[1].keys():
    data[key]=(prediction_stats[1][key] + prediction_stats[2][key])/2
  data["id"]=id
  return(data)

