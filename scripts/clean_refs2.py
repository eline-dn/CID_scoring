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
import json
from Bio.PDB import PDBParser, MMCIFParser, Superimposer, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, PDBIO, StructureBuilder
from Bio.PDB.Polypeptide import is_aa

from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer, MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import StructureBuilder
from Bio.PDB import PDBParser, Selection


### Path to this cloned GitHub repo:
SCRIPT_DIR = "/work/lpdi/users/eline/smol_binder_diffusion_pipeline"  # edit this to the GitHub repo path. Throws an error by default.
assert os.path.exists(SCRIPT_DIR)
sys.path.append(SCRIPT_DIR + "/scripts/utils")
import utils


#----------------------------------------------------------------------------------------------------------
"""-----------------------------------------------------SETUP-----------------------------------------------------"""
#----------------------------------------------------------------------------------------------------------

CONDAPATH = "/work/lpdi/users/eline/miniconda3"  # edit this depending on where your Conda environments live

### Path where the jobs will be run and outputs dumped
WDIR = "/work/lpdi/users/eline/smol_binder_diffusion_pipeline/1Z9Yout"
if not os.path.exists(WDIR):
    os.makedirs(WDIR, exist_ok=True)
print(f"Working directory: {WDIR}")
# Ligand information
LIGAND = "FUN"
MPNN_DIR = f"{WDIR}/1_proteinmpnn"

DIFFUSION_DIR = f"{WDIR}/0_diffusion"
AF3_DIR = f"{WDIR}/4_af3"
AF3_struct= f"{AF3_DIR}/output" ## run also with "output2"

### helper functions:
def copy_structure_with_only_chain(structure, chain_id):
  """
	From BindCraft's Biopyhton_utils : _copy_structure_with_only_chain (https://github.com/martinpacesa/BindCraft) 
	Return a new Structure containing only model 0 and a deep copy of chain `chain_id`."""
  # Build a tiny structure hierarchy: Structure -> Model(0) -> Chain(chain_id) -> Residues/Atoms

  sb = StructureBuilder.StructureBuilder()
  sb.init_structure("single")
  sb.init_model(1)
  sb.init_chain(chain_id)
  # Set segment ID, padded to 4 characters
  sb.init_seg(chain_id.ljust(4))  
  model0 = next(structure.get_models())
  if chain_id not in [c.id for c in model0.get_chains()]:
      raise ValueError(f"Chain '{chain_id}' not found.")
  chain = model0[chain_id]
  for res in chain:
      # Keep only amino-acid residues
      # Assuming is_aa is defined elsewhere and available
      if not is_aa(res, standard=False):
          continue
      hetflag, resseq, icode = res.id
      sb.init_residue(res.resname, hetflag, resseq, icode)

      for atom in res:
          sb.init_atom(atom.name, atom.coord, atom.bfactor, atom.occupancy,
                       atom.altloc, atom.fullname, element=atom.element)
  return sb.get_structure()

# change a chain's id:
def change_chain_id(structure,model_id, old_id, new_id, old_resname=None, new_resname=None):
  chain=structure[model_id][old_id]
  chain.id = new_id
  if new_resname is not None and old_resname is not None:
    for res in structure.get_residues():
      if res.resname.strip().startswith(old_resname):
        res.resname = new_resname
  return(structure)

def rm_ligands(pdb_str):
    return "\n".join(
        line for line in pdb_str.splitlines()
        if (line.startswith("ATOM"))
    ) + "\n"


MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

###-------

for confidence in glob.glob(f"{AF3_struct}/*/*_summary_confidences.json"):
    print("condidence:", confidence)
    design_name=os.path.basename(confidence)
    design_name=design_name.replace("_summary_confidences", "") # e.g. t2_2_7_4_t0.3_ltp0.2_dcut500.0
    design_name=design_name.replace(".json", "")

    # compute aligned rmsd to original rf diffusion bb:
    # retrieve bb: format t2_1_100_1_T0.3.pdb in MPNN_DIR
    bb_name = re.sub(r"_ltp0\.[1-9]_dcut(8\.0|15\.0|500\.0)", ".pdb", design_name)
    bb_name= bb_name.replace("t0","T0")
    bb_pdb=f"{MPNN_DIR}/backbones/{bb_name}"
    print("bb:", bb_pdb)
    # load / convert (?) model mmcif
    design_cif=confidence.replace("_summary_confidences.json","_model.cif")
    print("design_cif:", design_cif)
    
    pdb_parser = PDBParser(QUIET=True)
    mmcif_parser=MMCIFParser(QUIET= True)
    structure = pdb_parser.get_structure("ref", bb_pdb)
    print(bb_pdb)
    print(design_cif)
    
    model = structure[0]             # first model
    chain = model["A"]               # chain A
    # count only standard residues
    from Bio.PDB.Polypeptide import is_aa
    three_to_one = {
        "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I",
        "LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S",
        "THR":"T","VAL":"V","TRP":"W","TYR":"Y",
        # common variants
        "MSE":"M",  # Selenomethionine
    }
    residues = [res for res in chain.get_residues() if is_aa(res, standard=True)]
    #print("residues:", residues)
    print("len residues:", len(residues))
    # extract/trim the binder sequence
    binder_length=len(residues)-256
    print("binder len:", binder_length)
    
    from Bio.PDB import PDBIO, Chain
    from Bio.PDB.Polypeptide import is_aa
    # change chain id for binder residues (from A to B) and for ligand atoms: from B to L 
    # ligand:
    structure=change_chain_id(structure=structure,model_id=0, old_id="B", new_id="L", old_resname=None, new_resname=None)
    #binder:

    # Retrieve chain A (binder assumed to be first part)
    model = structure[0]             # first model
    chain_A = model["A"] 	    
    residues_A = list(chain_A.get_residues())
    #print("res A:", residues_A)
    if binder_length > len(residues_A):
        raise ValueError("binder_length exceeds number of residues in chain A")
    # Create new chain B
    chain_B = Chain.Chain("B")
    count=0
    # Transfer the binder residues into chain B
    for i, residue in enumerate(residues_A):
        if not is_aa(residue, standard=False):
            continue
        if i < binder_length:
            # Remove from chain A
            chain_A.detach_child(residue.id)
            # Add to chain B
            chain_B.add(residue)
            count+=1
    # Add new chain B to the model
    model.add(chain_B)
    print("len chain B:",count)
    io = PDBIO()
    io.set_structure(structure)
    io.save(f"/work/lpdi/users/eline/CID/1Z9Y_FUN/input/binder_refs/{design_name}.pdb")
