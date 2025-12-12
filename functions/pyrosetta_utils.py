"""
create ligand params from smiles and pdb for numbering (rdkit)
init with params
pdb to pose
cif to pose
atom contact count x 3
BC scoring (A-B)
BC scoring but with ligand

"""
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.MMCIFParser import *
import pyrosetta as pr
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.io import pose_from_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
import pyrosetta.distributed.io
import pyrosetta.rosetta.core.select.residue_selector as residue_selector
import os

import os
import math
import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree

# Minimal 3-letter -> 1-letter, including common alt codes
three_to_one_map = {
    "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I",
    "LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S",
    "THR":"T","VAL":"V","TRP":"W","TYR":"Y",
    # common variants
    "MSE":"M",  # Selenomethionine
}
# create ligand param file:
def create_param(pdb_file, smiles="c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl", lig_name='FUN'):
    pyrosetta.init(extra_options='-mute all') # required for test
    from rdkit_to_params import Params
    p = Params.from_smiles_w_pdbfile(pdb_file, smiles, name=lig_name) # the name has to match.
    print(p.is_aminoacid()) # False
    p.dump(f'{lig_name}.params')
    p.test().dump_pdb(f'parametrized_{lig_name}.pdb')
    print("Please ensure atoms names match between the param file and the pdb files you are going to analyse with Pyrosetta")
    return(f'{lig_name}.params')


# initialisation with a param file for ligand
# params = [f"FUN.params"] a list of string paths to param files

def pyr_init(params=None):
    # Rosetta params file(s)
    """
    Getting PyRosetta started
    """
    extra_res_fa = ""
    if params is not None:
        extra_res_fa = "-extra_res_fa"
        for p in params:
            extra_res_fa += f" {p}"
    
    NPROC = os.cpu_count()
    if "SLURM_CPUS_ON_NODE" in os.environ:
        NPROC = os.environ["SLURM_CPUS_ON_NODE"]
    elif "OMP_NUM_THREADS" in os.environ:
        NPROC = os.environ["OMP_NUM_THREADS"]
    
    
    DAB = f"{SCRIPT_PATH}/../utils/DAlphaBall.gcc" # This binary was compiled on UW systems. It may or may not work correctly on yours
    assert os.path.exists(DAB), "Please compile DAlphaBall.gcc and manually provide a path to it in this script under the variable `DAB`\n"\
                            "For more info on DAlphaBall, visit: https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/HolesFilter"
    
    
    pr.init(f"{extra_res_fa} -dalphaball {DAB} -beta_nov16 -run:preserve_header -mute all ")
            # f"-multithreading false -multithreading:total_threads {NPROC} -multithreading:interaction_graph_threads {NPROC}")
#---

def load_pose(pdbfile, assert_lig):
    pose = pyrosetta.pose_from_file(pdbfile)
    if assert_lig:
        ligand_resno = pose.size()
        if not pose.residue(ligand_resno).is_ligand():
            raise ValueError("the last residue is not a ligand")
    return(pose)

def load_pose_cif(cif_file, old_id="FUN", new_id="F", LIG_rename=None):
    from Bio.PDB import PDBIO
    import Bio.PDB
    CIF_parser = MMCIFParser(QUIET=True)
    structure = CIF_parser.get_structure("x", cif_file)
    chain=structure[0][old_id]
    chain.id = new_id
    if LIG_rename is not None:
        for res in structure.get_residues():
            if res.resname.strip().startswith("LIG"):  # or your rule
                res.resname = LIG_rename
    io = PDBIO()
    io.set_structure(structure)
    pdbfile=cif_file.replace(".cif", ".pdb") # convert mmcif to pdb
    io.save(pdbfile)
    pose = pyrosetta.pose_from_file(pdbfile)
    os.remove(pdbfile)  # tidy up pdb
    return(pose)


# compute contacts between two chains based on Van der Waals radii (rVdW + 0.2 Å tolerance) of each pair of atoms
"""  draft
def contact_counts(pdbfile, ch_id1, ch_id2):
    pose = pyrosetta.pose_from_file(pdbfile)
    scorefxn = pyrosetta.get_fa_scorefxn()
    totE=scorefxn(pose)
    graph = pose.energies().tenA_neighbor_graph()
  
#pose.residue(i).atom(j).atom_type().lj_radius()
#distance <= lj_radius(i) + lj_radius(j) + 0.2
In PyRosetta you retrieve LJ radii from:

atom_type_set = pose.residue(i).atom_type_set()
atom_type = pose.residue(i).atom_type(atom_index)
lj_radius = atom_type.lj_radius()
"""

# one option:

import math

def count_chain_atom_contacts(pose, ch_id1, ch_id2, delta=0.2):

    # ensure neighbor graph is built
    sf = pyrosetta.get_fa_scorefxn()
    sf(pose)

    graph = pose.energies().tenA_neighbor_graph()

    total_contacts = 0

    edge_it = graph.const_edge_list_begin()
    edge_end = graph.const_edge_list_end()

    while edge_it != edge_end:
        edge = edge_it.__deref__()

        i = edge.lower_node_ind()
        j = edge.upper_node_ind()

        chain_i = pose.pdb_info().chain(i)
        chain_j = pose.pdb_info().chain(j)

        # only consider residues from the two target chains
        if {chain_i, chain_j} != {ch_id1, ch_id2}:
            edge_it.increment()
            continue

        res_i = pose.residue(i)
        res_j = pose.residue(j)

        # iterate atoms of the two residues
        for ai in range(1, res_i.natoms() + 1):
            xyz_i = res_i.xyz(ai)
            lj_i = res_i.atom_type(ai).lj_radius()

            for aj in range(1, res_j.natoms() + 1):
                xyz_j = res_j.xyz(aj)
                lj_j = res_j.atom_type(aj).lj_radius()

                # squared distance
                dsq = xyz_i.distance_squared(xyz_j)
                cutoff = (lj_i + lj_j + delta)

                # compare using squared values
                if dsq < cutoff * cutoff:
                    total_contacts += 1

        edge_it.increment()

    return total_contacts
    
# another option
def get_atomic_contact_data(pose):
    '''Get the atomic contact number for a pose.'''
    acf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<AtomicContactCount name="contact" />')
    acf.initialize_cross_chain()
    return [(structure_name, acf.report_sm(pose))]

# third option:
def count_atomic_contacts(pdb, chain1, chain2, delta=0.2):
    VdWradii={"N":1.55, "C":1.70, "O":1.52, "H":1.20,"S":1.80, "CL": 1.75}
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb)
    
    if (chain1 not in [c.id for c in structure[0].get_chains()]) or (chain2 not in [c.id for c in structure[0].get_chains()]):
            raise ValueError(f"Chain not found.")
    # Get the specified chain
    binder_atoms = Selection.unfold_entities(structure[0][chain1], 'A') # get all the child atoms from the residues in chain 1 
    binder_coords = np.array([atom.coord for atom in binder_atoms])
    target_atoms = Selection.unfold_entities(structure[0][chain2], 'A') # get all the child atoms from the residues in chain 2
    target_coords = np.array([atom.coord for atom in target_atoms])
   
    # Build KD trees for both chains
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)
    
    atom_distance_cutoff=1.6 *2 + delta #appreoxiamtion to reproduce neosurf distance
    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)
    return(len(pairs))

# find hotspot residues:
# identify interacting residues at the binder interface
# from bindcraft, modified in order to be able to add the ligand in the equation

def hotspot_residues(trajectory_pdb, binder_chain="B", target_plus_lig=False, lig_only=False, lig_chain="F", atom_distance_cutoff=4.0):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", trajectory_pdb)
    if target_plus_lig==True and lig_only==True:
        raise ValueError(f"do not use target_plus_lig and lig_only at the same time!")
    if target_plus_lig==True or lig_only==True:
        if lig_chain not in [c.id for c in structure[0].get_chains()]:
            raise ValueError(f"Chain '{lig_chain}' not found. Check ligand chain!")
    # Get the specified chain
    binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A') # get all the child atoms from the residues in  binder  
    binder_coords = np.array([atom.coord for atom in binder_atoms])
    target_atoms = Selection.unfold_entities(structure[0]['A'], 'A') # get all the child atoms from the residues in chain A
    ligand_atoms=Selection.unfold_entities(structure[0][lig_chain], 'A') # get all the child atoms from the residues in chain A
    
    if target_plus_lig:       
        target_atoms+=ligand_atoms
        target_coords = np.array([atom.coord for atom in target_atoms]) # atom coordinate for target protein + ligand
    elif lig_only:
        target_coords = np.array([atom.coord for atom in ligand_atoms]) # atom coordinate for ligand only
    else: # atom coordinate for target protein only
        target_coords = np.array([atom.coord for atom in target_atoms])
    # Build KD trees for both chains
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    # Prepare to collect interacting residues
    interacting_residues = {}

    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    # Process each binder atom's interactions
    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        # Convert three-letter code to single-letter code using the manual dictionary
        if binder_resname in three_to_one_map:
            aa_single_letter = three_to_one_map[binder_resname]
            for close_idx in close_indices:
                target_residue = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter
    return interacting_residues

# covert cif to pdb
def cif2pdb(cif_path, pdb_path):
  parser = MMCIFParser(QUIET=True)
  structure = parser.get_structure("x", cif_path)
  io = PDBIO(use_model_flag=1)
  io.set_structure(structure)
  io.save(pdb_path)


# BindCraft's scoring function from pdb


# Rosetta interface scores
def score_nolig_interface(pdb_file, interface= "A_B", binder_chain="B"):
    # load pose
    pose = pr.pose_from_pdb(pdb_file)
    # analyze interface statistics
    iam = InterfaceAnalyzerMover()
    iam.set_interface(interface)
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

    # Initialize dictionary with all amino acids
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}

    # Initialize list to store PDB residue IDs at the interface
    interface_residues_set = hotspot_residues(pdb_file, binder_chain, target_plus_lig=False, lig_only=False, lig_chain="F")
    interface_residues_pdb_ids = []
    
    # Iterate over the interface residues
    for pdb_res_num, aa_type in interface_residues_set.items():
        # Increase the count for this amino acid type
        interface_AA[aa_type] += 1

        # Append the binder_chain and the PDB residue number to the list
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

    # count interface residues
    interface_nres = len(interface_residues_pdb_ids)

    # Convert the list into a comma-separated string
    interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)

    # Calculate the percentage of hydrophobic residues at the interface of the binder
    hydrophobic_aa = set('ACFILMPVWY')
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    if interface_nres != 0:
        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
    else:
        interface_hydrophobicity = 0

    # retrieve statistics
    interfacescore = iam.get_all_data()
    interface_sc = interfacescore.sc_value # shape complementarity
    interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
    interface_dG = iam.get_interface_dG() # interface dG
    interface_dSASA = iam.get_interface_delta_sasa() # interface dSASA (interface surface area)
    interface_packstat = iam.get_interface_packstat() # interface pack stat score
    interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100 # ratio of dG/dSASA (normalised energy for interface area size)
    buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100 # Hbonds per interface size percentage
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100 # Unsaturated H-bonds per percentage
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None

    # calculate binder energy score
    chain_design = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)

    # calculate binder SASA fraction
    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)

    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    else:
        interface_binder_fraction = 0

    # calculate surface hydrophobicity
    binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]

    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core = False, pick_boundary = False, pick_surface = True)
    surface_res = layer_sel.apply(binder_pose)

    exp_apol_count = 0
    total_count = 0 
    
    # count apolar and aromatic residues at the surface
    for i in range(1, len(surface_res) + 1):
        if surface_res[i] == True:
            res = binder_pose.residue(i)

            # count apolar and aromatic residues as hydrophobic
            if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
                exp_apol_count += 1
            total_count += 1

    surface_hydrophobicity = exp_apol_count/total_count

    # output interface score array and amino acid counts at the interface
    interface_scores = {
    'binder_score': binder_score,
    'surface_hydrophobicity': surface_hydrophobicity,
    'interface_sc': interface_sc,
    'interface_packstat': interface_packstat,
    'interface_dG': interface_dG,
    'interface_dSASA': interface_dSASA,
    'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
    'interface_fraction': interface_binder_fraction,
    'interface_hydrophobicity': interface_hydrophobicity,
    'interface_nres': interface_nres,
    'interface_interface_hbonds': interface_interface_hbonds,
    'interface_hbond_percentage': interface_hbond_percentage,
    'interface_delta_unsat_hbonds': interface_delta_unsat_hbonds,
    'interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage
    }

    # round to two decimal places
    interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

    return interface_scores, interface_AA, interface_residues_pdb_ids_str



# try to score the binder vs ligand + target interface:
# Rosetta interface scores
def score_lig_interface(pdb_file, interface= "AF_B", binder_chain="B"):
    # load pose
    pose = pr.pose_from_pdb(pdb_file)
    # analyze interface statistics
    iam = InterfaceAnalyzerMover()
    iam.set_interface(interface)
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

    # Initialize dictionary with all amino acids
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}

    # Initialize list to store PDB residue IDs at the interface
    interface_residues_set = hotspot_residues(pdb_file, binder_chain, target_plus_lig=True, lig_only=False, lig_chain="F")
    interface_residues_pdb_ids = []
    interface_residues_set_lig_only = hotspot_residues(pdb_file, binder_chain, target_plus_lig=False, lig_only=True, lig_chain="F")
    
    # Iterate over the interface residues
    for pdb_res_num, aa_type in interface_residues_set.items():
        # Increase the count for this amino acid type
        interface_AA[aa_type] += 1

        # Append the binder_chain and the PDB residue number to the list
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

    # count interface residues
    interface_nres = len(interface_residues_pdb_ids)
    interface_nres_lig_only=len(interface_residues_set_lig_only)
    
    # Convert the list into a comma-separated string
    interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)

    # Calculate the percentage of hydrophobic residues at the interface of the binder
    hydrophobic_aa = set('ACFILMPVWY')
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    if interface_nres != 0:
        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
        interface_lig_prop=(interface_nres_lig_only/interface_nres) *100
    else:
        interface_hydrophobicity = 0
        interface_lig_prop=0

    # retrieve statistics
    interfacescore = iam.get_all_data()
    interface_sc = interfacescore.sc_value # shape complementarity
    interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
    interface_dG = iam.get_interface_dG() # interface dG
    interface_dSASA = iam.get_interface_delta_sasa() # interface dSASA (interface surface area)
    interface_packstat = iam.get_interface_packstat() # interface pack stat score
    interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100 # ratio of dG/dSASA (normalised energy for interface area size)
    buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100 # Hbonds per interface size percentage
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100 # Unsaturated H-bonds per percentage
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None

    # calculate binder energy score
    chain_design = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)

    # calculate binder SASA fraction
    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)

    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    else:
        interface_binder_fraction = 0

    # calculate surface hydrophobicity
    binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]

    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core = False, pick_boundary = False, pick_surface = True)
    surface_res = layer_sel.apply(binder_pose)

    exp_apol_count = 0
    total_count = 0 
    
    # count apolar and aromatic residues at the surface
    for i in range(1, len(surface_res) + 1):
        if surface_res[i] == True:
            res = binder_pose.residue(i)

            # count apolar and aromatic residues as hydrophobic
            if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
                exp_apol_count += 1
            total_count += 1

    surface_hydrophobicity = exp_apol_count/total_count

    # output interface score array and amino acid counts at the interface
    interface_scores = {
    'wlig_binder_score': binder_score,
    'wlig_surface_hydrophobicity': surface_hydrophobicity,
    'wlig_interface_sc': interface_sc,
    'wlig_interface_packstat': interface_packstat,
    'wlig_interface_dG': interface_dG,
    'wlig_interface_dSASA': interface_dSASA,
    'wlig_interface_dG_SASA_ratio': interface_dG_SASA_ratio,
    'wlig_interface_fraction': interface_binder_fraction,
    'wlig_interface_hydrophobicity': interface_hydrophobicity,
    'wlig_interface_nres': interface_nres,
    'binder_lig_interface_nres':interface_nres_lig_only,
    'binder_lig_interface_%':interface_lig_prop,
    'wlig_interface_interface_hbonds': interface_interface_hbonds,
    'wlig_interface_hbond_percentage': interface_hbond_percentage,
    'wlig_interface_delta_unsat_hbonds': interface_delta_unsat_hbonds,
    'wlig_interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage
    }

    # round to two decimal places
    interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

    return interface_scores, interface_AA, interface_residues_pdb_ids_str
