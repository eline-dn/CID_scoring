"""
biopython utils:
load from cif
load from pdb
extract sequence from chain
extract chain length
align to chain and RMSD
RMSD no alignment 
write pdb
write cif
extract chain from struct -> struct
"""

### Import dependencies
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

# Minimal 3-letter -> 1-letter, including common alt codes
three_to_one_map = {
    "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I",
    "LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S",
    "THR":"T","VAL":"V","TRP":"W","TYR":"Y",
    # common variants
    "MSE":"M",  # Selenomethionine
}

# load from PDB:
def load_PDB(pdb_file):
  parser = PDBParser(QUIET=True)
  structure = parser.get_structure("x", pdb_file)
  return(structure)

# load from CIF:

def load_CIF(cif_file):
  parser = MMCIFParser(QUIET=True)
  structure = parser.get_structure("x", cif_file)
  return(structure)

# extract sequence from chain 
def chain2seq(structure,  chain_id:str, make_str=True):
  chain=structure[0][chain_id]
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

# extract chain length
def chain2length(structure, chain_id:str):
  chain=structure[0][chain_id]
  # count only standard residues
  residues = [res for res in chain.get_residues() if res.id[0] == " "] #res.id[0] == " " filters out HETATM residues.
  return(len(residues))

# change a chain's id:
def change_chain_id(structure,model_id, old_id, new_id, old_resname=None, new_resname=None):
  chain=structure[model_id][old_id]
  chain.id = new_id
  if new_resname is not None and old_resname is not None:
    for res in structure.get_residues():
      if res.resname.strip().startswith(old_resname):
        res.resname = new_resname
  return(structure)

# focus on aligning only one chain from the pdbs (e.g. a target chain) and then return this chain's rmsd or the full structure's rmsd

# Helper: get ordered CA coordinate lists based on sequence
def get_ca_atoms(chain):
    atoms = []
    for res in chain:
        if is_aa(res, standard=True) and "CA" in res:
            atoms.append(res["CA"])
    return atoms

# Helper: compute aligned rmsd with superimposer from ca list obtained with get_ca_atoms
def sup_rmsd(ref_ca, mov_ca): # compute aligned rmsd with superimposer from ca list obtained with get_ca_atoms
  # ensure same length by trimming the longer one
  L = min(len(ref_ca), len(mov_ca))
  ref_ca = ref_ca[:L]
  mov_ca = mov_ca[:L]

  sup = Superimposer()
  sup.set_atoms(ref_ca, mov_ca)
  rot, tran = sup.rotran
  rmsd = sup.rms
  return(rmsd)
  
## /!\ only for amino acid chains!   
def aligned_chain_rmsd(ref, mov, ref_chain_id, mov_chain_id): # compute one aligned chain rmsd from structures
  ref_chain=ref[0][ref_chain_id]
  mov_chain=mov[0][mov_chain_id]
  ref_ca = get_ca_atoms(ref_chain)
  mov_ca = get_ca_atoms(mov_chain)
  rmsd=sup_rmsd(ref_ca, mov_ca)
  return rmsd
 
def full_aligned_rmsd(ref, mov, mapping):
  # concat all the chains' CA in the same order
  # provide in this case a mapping of the ref and mov structure chains in the shape
  # { ref_chain_id1: mov_chain_id1, ref_chain_id2: mov_chain_id2,...}
  ref_ca_list=list()
  mov_ca_list=list()
  for ref_chain_id, mov_chain_id in mapping.items():
    ref_chain=ref[0][ref_chain_id]
    mov_chain=mov[0][mov_chain_id]
    ref_ca = get_ca_atoms(ref_chain)
    mov_ca = get_ca_atoms(mov_chain)
    ref_ca_list+=ref_ca
    mov_ca_list+=mov_ca
  rmsd=sup_rmsd(ref_ca_list, mov_ca_list)
  return rmsd

def unaligned_rmsd(ref, mov, mapping): # compute unaligned rmsd for the chains in mapping, not for ligands!
  def ca_map(chain):
    out = {}
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        if "CA" in res:
            hetflag, resseq, icode = res.get_id()
            out[(resseq, icode)] = res["CA"]
    return out
  
  ref_ca_list={}
  mov_ca_list={}
  for ref_chain_id, mov_chain_id in mapping.items():
    ref_chain=ref[0][ref_chain_id]
    mov_chain=mov[0][mov_chain_id]
    ref_ca = ca_map(ref_chain)
    mov_ca = ca_map(mov_chain)
    ref_ca_list={**ref_ca_list, **ref_ca}
    mov_ca_list={**mov_ca_list, **mov_ca}
  
  # Common residue keys, ordered by residue number then insertion code
  common = sorted(set(ref_ca_list.keys()).intersection(mov_ca_list.keys()),
                  key=lambda k: (k[0], (k[1] or " ")))
  
  if len(common) < 3:
      raise ValueError(
          f"Not enough matched residues with Cα to compute RMSD without alignment "
          f"(found {len(common)})."
      )
  
  ref_coords = np.array([ref_ca_list[k].get_coord() for k in common], dtype=float)
  mov_coords = np.array([mov_ca_list[k].get_coord() for k in common], dtype=float)
  
  # Unaligned RMSD = sqrt(mean(||ref - mov||^2))
  diffs = ref_coords - mov_coords
  rmsd = float(np.sqrt((diffs * diffs).sum(axis=1).mean()))
  print(round(rmsd, 2))
  return(rmsd)

def unaligned_ligand_rmsd(ref, mov,ref_lig_chain_id,mov_lig_chain_id):
  # ligand atoms: compute RMSD over all heavy atoms
  ref_ligand=ref[0][ref_lig_chain_id]
  mov_ligand=ref[0][mov_lig_chain_id]
  ref_lig_atoms = np.array([a.get_coord() for a in ref_ligand.get_atoms()])
  mov_lig_atoms = np.array([a.get_coord() for a in mov_ligand.get_atoms()])
  L = min(len(ref_lig_atoms), len(mov_lig_atoms))
  diffs=ref_lig_atoms[:L]-mov_lig_atoms[:L]
  lig_rmsd = float(np.sqrt((diffs * diffs).sum(axis=1).mean()))
  return(lig_rmsd)

def align_to_chain(ref,mov,mapping): # align a structure to a reference structure based on the chains in mapping
  # concat all the chains' CA in the same order
  # provide in this case a mapping of the ref and mov structure chains in the shape
  # { ref_chain_id1: mov_chain_id1, ref_chain_id2: mov_chain_id2,...}
  # !!! the targeted chains should be the same lentgh!!
  ref_ca_list=list()
  mov_ca_list=list()
  for ref_chain_id, mov_chain_id in mapping.items():
    ref_chain=ref[0][ref_chain_id]
    mov_chain=mov[0][mov_chain_id]
    ref_ca = get_ca_atoms(ref_chain)
    mov_ca = get_ca_atoms(mov_chain)
    ref_ca_list+=ref_ca
    mov_ca_list+=mov_ca

  if len(ref_ca_list) != len(mov_ca_list):
    raise ValueError(
        f" the list of atoms to align are of different lengths, check the chain mapping and your chain lengths"
        f"(number of ref atoms {len(ref_ca_list)} and number of mov atoms {len(mov_ca_list)} )."
    )
  sup = Superimposer()
  sup.set_atoms(ref_ca_list, mov_ca_list)
  rot, tran = sup.rotran
  # Apply transform to ALL atoms in the moving structure
  for atom in mov.get_atoms():
      atom.transform(rot, tran)
  # return the modified mov structure: (the structure will be modified anyway, even if you use a copy after that
  return(mov)


def write_pdb(structure, output_pdb_path):
  io = PDBIO(use_model_flag=1)
  io.set_structure(structure)
  io.save(output_pdb_path)

def write_cif(structure, output_cif_path):
  io=MMCIFIO()
  io.set_structure(structure)
  io.save(output_cif_path)

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
  model0 = structure[0]
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

#### some higher level wrapper functions to prepare the reference pdbs and ensure a coherent chain notation and numbering
# chain names:
# numbering
# checl ligand notation

#### some higher level wrapper functions to analyse the structure repredictions and comparing them to the references
## in the case of a ternary complex:
def ternary_RMSDs(ref, mov, mapping):
	data={}
	data["target_RMSD"]=aligned_chain_rmsd(ref, mov, "A", "A")
	data["binder_RMSD"]=aligned_chain_rmsd(ref, mov, "B", "B")
	data["full_ RMSD"]=full_aligned_rmsd(ref, mov, mapping)
	mov_aligned=align_to_chain(ref,mov,mapping)
	data["d_binding_site"]=unaligned_rmsd(ref, mov_aligned, {"B":"B"})
	ref_lig_chain_id=list(mapping.keys())[-1]
	mov_lig_chain_id=mapping[ref_lig_chain_id]
	data["ligand_rmsd"]=unaligned_ligand_rmsd(ref, mov_aligned, ref_lig_chain_id ,mov_lig_chain_id)
	return(data)

def binary_RMSDs(ref, mov, mapping):
	data={}
	data["target_RMSD"]=aligned_chain_rmsd(ref, mov, "A", "A")
	data["binder_RMSD"]=aligned_chain_rmsd(ref, mov, "B", "B")
	data["full_ RMSD"]=full_aligned_rmsd(ref, mov, mapping)
	mov_aligned=align_to_chain(ref,mov,mapping)
	data["d_binding_site"]=unaligned_rmsd(ref, mov_aligned, {"B":"B"})
	return(data)

def af3_out_2_norm_pdb(cif_path, lig=True, lig_name=None):
	structure=load_CIF(cif_path)
	pdb_path=cif_path.replace(".cif", ".pdb")
	
	first_model = next(structure.get_models())
    first_model_id = first_model.id
    
	if lig: # if a ligand is present
		if lig_name is None:
			raise ValueError("Specify lig_name argument if the cif file contains a ligand")
		change_chain_id(structure,first_model_id, lig_name, "L", "LIG", lig_name) # replace chain id "FUN" with "L" and ligand's resname with its name and not e.g. LIG_FUN

	write_pdb(structure, pdb_path)
	return(pdb_path)

