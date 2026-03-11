from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import os, shutil
import pandas as pd
import numpy as np
import os
import math
import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer, MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import StructureBuilder
from Bio.PDB import PDBParser, Selection
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
import numpy as np

three_to_one_map = {
    "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I",
    "LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S",
    "THR":"T","VAL":"V","TRP":"W","TYR":"Y",
    # common variants
    "MSE":"M",  # Selenomethionine
}

one_to_three_map = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR"
}


def compute_side_chain_plddt(pdb_path, residue_number, residue_code):
    """
    Compute the mean side-chain pLDDT for a given residue in chain B of the first model.

    Args:
        pdb_path (str): Path to the PDB file.
        residue_number (int): Residue number (e.g., 123).
        residue_code (str): Three-letter residue code (e.g., "ALA").

    Returns:
        float: Mean side-chain pLDDT, or None if residue not found or no side-chain atoms.
    """
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_path)
    # Get the first model and chain B
    model = next(structure.get_models())
    chain = model["B"]
    # Find the residue
    residue = None
    for res in chain:
        print(res.get_id()[1], res.get_resname())
        if res.get_id()[1] == residue_number and res.get_resname() == residue_code:
            print("found")
            residue = res
            break
    if residue is None:
        print(f"Residue {residue_code} {residue_number} not found in chain B.")
        return None
    # Collect side-chain atoms (not CA, N, O, C)
    side_chain_atoms = []
    for atom in residue:
        if atom.get_name() not in ("CA", "N", "O", "C"):
            side_chain_atoms.append(atom)
    if not side_chain_atoms:
        print(f"No side-chain atoms found for residue {residue_code} {residue_number}.")
        return None
    # Sum pLDDT values
    plddt_sum = 0.0
    for atom in side_chain_atoms:
        # pLDDT is stored in the B-factor field (atom.get_bfactor())
        plddt_sum += atom.get_bfactor()
    mean_plddt = plddt_sum / len(side_chain_atoms)
    return mean_plddt



# Example usage:
# mean_plddt = compute_side_chain_plddt("example.pdb", 123, "ALA")
# print(f"Mean side-chain pLDDT: {mean_plddt}")
def load_PDB(pdb_file):
  parser = PDBParser(QUIET=True)
  structure = parser.get_structure("x", pdb_file)
  return(structure)

def chain2seq(structure,  chain_id:str, make_str=True):
  chain=next(structure.get_models())[chain_id]
  from Bio.PDB.Polypeptide import is_aa
  # count only standard residues
  residues = [res for res in chain.get_residues() if is_aa(res, standard=True)]
  # convert residues to one-letter sequence
  res_letters = []
  for res in residues:
      try:
          aa = three_to_one_map[res.resname]
          res_letters.append(aa)
      except KeyError:
          raise ValueError(f"Unknown residue: {res.resname} at {res.id}")
  if make_str==False:
    return res_letters
  else:
    return("".join(res_letters))


import glob
def find_pdb(id, outdir):
  pdb_path=os.path.join(outdir +"/"+ id, id+"_model.pdb")
  #pdb_path=glob.blop
  if os.path.isfile(pdb_path):
    return pdb_path
  else:
    print(Warning(f"PDB file not found: {pdb_path} for id {id}"))
    return None

df=pd.read_csv("/work/lpdi/users/eline/CID/b2_FUN_1Z9Y/output/pass_af3/pass_pyr/cid-hbonding-aa.csv")

#outdir="/work/lpdi/users/eline/CID/b2_FUN_1Z9Y/output/pass_af3/pass_pyr" # batch 2
#outdir="/work/lpdi/users/eline/CID/FUN_1Z97/output/af3ternary/accepted3" # batch 1
#df["sc_plddt"]=np.nan
#df["binder_sequence"]=""
for ID in df.Name:
  """aa = df.loc[df['Name'] == ID, "H_bonding_aa"].values[0]
  res_num=int(aa[:-1])
  res_code=str(aa[-1])"""
  path=find_pdb(ID, outdir)
  if path is not None:
    #sc_plddt=compute_side_chain_plddt(path, res_num, one_to_three_map[res_code])
    struct=load_PDB(path)
    binder_sequence=chain2seq(struct, chain_id="B", make_str=True)
    df.loc[df['Name'] == ID, "binder_sequence"]=binder_sequence
  else:
    print(f"skipping id {ID}, file not found")

df.to_csv("/work/lpdi/users/eline/CID/b2_FUN_1Z9Y/output/pass_af3/pass_pyr/cid_designs.csv", index=False)

import matplotlib.pyplot as plt
"""
# Example: Plot histogram for column 'sc_plddt' in DataFrame df
plt.figure(figsize=(8, 5))
plt.hist(df['sc_plddt'], bins=20, edgecolor='black', alpha=0.7)
plt.title('Histogram of sc_plddt Values')
plt.xlabel('sc_plddt')
plt.ylabel('Frequency')
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()"""