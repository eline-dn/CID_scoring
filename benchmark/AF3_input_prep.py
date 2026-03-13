import os
import csv
import json
import sys
import copy
import shutil
from Bio.PDB import MMCIFIO, Select
from Bio.PDB import PDBList
# import AF3_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from biopython_utils import load_CIF, chain2seq, change_chain_id, copy_structure_with_only_chain

# --- Configuration ---
CSV_INPUT = "data_complexes.csv"  # columns: pdb_id, chain_a, chain_b, ligand_id, smiles
OUTPUT_BASE_DIR = "af3_inputs"
PDB_STORAGE_DIR = "pdb_files"  # Where original pdb/cif files are stored/downloaded

""" usage:
AF3_input_prep
"""
CIF_HEADER_CONTENT = """#
loop_
_pdbx_audit_revision_history.ordinal 
_pdbx_audit_revision_history.data_content_type 
_pdbx_audit_revision_history.major_revision 
_pdbx_audit_revision_history.minor_revision 
_pdbx_audit_revision_history.revision_date 
1 'Structure model' 1 0 2006-05-23 
2 'Structure model' 1 1 2008-04-30 
3 'Structure model' 1 2 2011-07-13 
4 'Structure model' 1 3 2017-10-11 
# 
"""

three_to_one = {
    "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I",
    "LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S",
    "THR":"T","VAL":"V","TRP":"W","TYR":"Y",
    # common variants
    "MSE":"M",  # Selenomethionine
}

class ChainSelect(Select):
    """Filter to save only a specific chain."""
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

def ensure_cif_exists(pdb_id, storage_dir):
    """
    Checks if a CIF file exists locally; if not, downloads it from RCSB.
    Returns the path to the file.
    """
    if not os.path.exists(storage_dir):
        os.makedirs(storage_dir)
    
    # Biopython downloads files as 'pdbXXXX.ent' or 'XXXX.cif' 
    # and often uses lowercase for the filename.
    expected_file = os.path.join(storage_dir, f"{pdb_id.lower()}.cif")
    downloaded_file = expected_file  # Initialize to expected path
    
    if not os.path.exists(expected_file):
        print(f"File {expected_file} not found. Downloading from RCSB...")
        pdbl = PDBList()
        # file_format='mmCif' ensures we get the .cif version
        pdbl.retrieve_pdb_file(pdb_id, pdir=storage_dir, file_format='mmCif')
        
        # PDBList often names the file 'pdbID.cif', let's standardize it
        if not os.path.exists(downloaded_file):
            # Handle cases where PDBList might name it 'pdbXXXX.cif'
            alt_name = os.path.join(storage_dir, f"pdb{pdb_id.lower()}.cif")
            if os.path.exists(alt_name):
                os.rename(alt_name, downloaded_file)
    
    return downloaded_file

def save_cif_with_header(structure, chain_id, output_path):
    """Saves a single chain to CIF and injects the AF3 header after the data_ line.

    AF3 requires the CIF to start with a data_ record.  The Biopython save
    produces a file where the first line is already the data_ header, so
    we must insert our custom CIF_HEADER_CONTENT *after* that line instead of
    prepending it.  Prepending would push the data_ line down and trigger
    "INVALID_ARGUMENT: The CIF file does not start with the data_ field."
    when AF3 loads the template.
    """
    io = MMCIFIO()
    io.set_structure(structure)
    temp_path = output_path + ".tmp"
    io.save(temp_path, ChainSelect(chain_id))

    # read the saved file and insert header after first data_ line
    with open(temp_path, 'r') as original:
        lines = original.readlines()
    if not lines:
        raise ValueError(f"Empty CIF produced for chain {chain_id}")
    # find the first line starting with data_
    if lines[0].startswith("data_"):
        # keep the first line, then add our header, then the rest
        new_lines = [lines[0]] + [CIF_HEADER_CONTENT] + lines[1:]
    else:
        # fallback: just prepend (shouldn't usually happen)
        new_lines = [CIF_HEADER_CONTENT] + lines

    with open(output_path, 'w') as final:
        final.writelines(new_lines)
    os.remove(temp_path)

def create_af3_json(name, seq_a, path_a, seq_b, path_b, ligand_smiles=None, ligand_id=None):
    """Constructs the AF3 JSON object with templates for both chains."""
    indices_a = list(range(len(seq_a)))
    indices_b = list(range(len(seq_b)))
    
    # Protein A (with Template)
    data = {
        "name": name,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": seq_a,
                    "unpairedMsa": "",
                    "pairedMsa": "",
                    "templates": [
                        {
                            "mmcifPath": os.path.abspath(path_a),
                            "queryIndices": indices_a,
                            "templateIndices": indices_a
                        }
                    ]
                }
            },
            {
                "protein": {
                    "id": ["B"],
                    "sequence": seq_b,
                    "unpairedMsa": "",
                    "pairedMsa": "",
                    "templates": [
                        {
                            "mmcifPath": os.path.abspath(path_b),
                            "queryIndices": indices_b,
                            "templateIndices": indices_b
                        }
                    ]
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }

    # Optional Ligand (Ternary Complex)
    if ligand_smiles and ligand_id:
        data["sequences"].append({
            "ligand": {
                "id": ligand_id,
                "smiles": ligand_smiles
            }
        })
    
    return data

def sanitize_ligand_id(ligand_id: str) -> str:
    """Return a single uppercase-letter ligand ID safe for AF3.

    AF3 requires ligand identifiers to be uppercase letters.  If the provided
    ID is empty or not a single uppercase letter, fall back to 'L'.  If the
    CSV already contains numeric or multi-character IDs, we skip numeric
    characters and take the first alphabetic character, uppercased. If none
    found, default to 'L'.
    """
    if ligand_id is None:
        return "L"
    ligand = ligand_id.strip()
    if ligand.isalpha() and ligand.isupper():
        return ligand
    # No alphabetic character found, use default
    return "L"


def count_chains(structure):
    """Count the number of chains in a structure."""
    model = next(structure.get_models())
    return len(list(model.get_chains()))

def has_hetatm_records(structure):
    """Check if a structure contains any HETATM records (non-standard residues)."""
    for residue in structure.get_residues():
        # res.id[0] != " " indicates a HETATM record
        if residue.id[0] != " ":
            return True
    return False

def validate_template(structure, expected_chain_id, template_name):
    """
    Validate a template structure.
    Returns (is_valid, error_message)
    
    Checks:
    - Structure contains exactly one chain
    - Structure contains no HETATM records
    - Chain ID matches expected ID (A or B)
    """
    # Check number of chains
    num_chains = count_chains(structure)
    if num_chains != 1:
        return False, f"{template_name} has {num_chains} chains, expected 1"
    
    # Check for HETATM records
    if has_hetatm_records(structure):
        return False, f"{template_name} contains HETATM records"
    
    # Check chain ID
    model = next(structure.get_models())
    chain = next(model.get_chains())
    if chain.get_id() != expected_chain_id:
        return False, f"{template_name} has chain '{chain.get_id()}', expected '{expected_chain_id}'"
    
    return True, None

def validate_starting_pdb(structure, chain_id, pdb_id):
    """
    Validate the starting PDB for a specified chain.
    Returns (is_valid, error_message)
    
    Checks:
    - Chain exists in structure
    - Chain has at least one standard amino acid after filtering HETATM
    """
    model = next(structure.get_models())
    
    # Check if chain exists
    try:
        chain = model[chain_id]
    except KeyError:
        return False, f"Chain '{chain_id}' not found in {pdb_id}"
    
    # Check that there is at least one standard amino acid after filtering HETATM
    standard_residues = [res for res in chain.get_residues() if res.id[0] == " "]
    if not standard_residues:
        return False, f"Chain '{chain_id}' in {pdb_id} has no standard amino acids after filtering HETATM"
    
    return True, None

def remove_hetatm_from_structure(structure):
    """Remove HETATM records (non-standard residues) from a structure."""
    model = next(structure.get_models())
    residues_to_remove = [res for res in model.get_residues() if res.id[0] != " "]
    for residue in residues_to_remove:
        chain = residue.get_parent()
        chain.detach_child(residue.id)

def process_complexes(csv_file):
    if not os.path.exists(OUTPUT_BASE_DIR):
        os.makedirs(OUTPUT_BASE_DIR)
    
    # Create directories for validated JSONs
    json_bin_dir = os.path.join(OUTPUT_BASE_DIR, "json_bin")
    json_ter_dir = os.path.join(OUTPUT_BASE_DIR, "json_ter")
    os.makedirs(json_bin_dir, exist_ok=True)
    os.makedirs(json_ter_dir, exist_ok=True)

    with open(csv_file, mode='r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_id = row['pdb_id']
            chain_a_id = row['chain_a']
            chain_b_id = row['chain_b']
            ligand_id = sanitize_ligand_id(row['ligand_id'])
            smiles = row['smiles']

            # Skip rows with empty smiles
            if not smiles or smiles.strip() == '':
                print(f"Skipping {pdb_id}: empty smiles field.")
                continue

            # Setup directory for this complex
            work_dir = os.path.join(OUTPUT_BASE_DIR, f"{pdb_id}_{chain_a_id}_{chain_b_id}")
            os.makedirs(work_dir, exist_ok=True)
            
            # retrieve CIF file (download if not present)
            try:
                cif_file = ensure_cif_exists(pdb_id, PDB_STORAGE_DIR)
            except Exception as e:
                print(f"Error retrieving {pdb_id}: {e}, skipping.")
                continue
            
            # Load Structure
            structure = load_CIF(cif_file)
            
            # Validate starting PDB for both chains
            is_valid_a, error_a = validate_starting_pdb(structure, chain_a_id, pdb_id)
            if not is_valid_a:
                print(f"Skipping {pdb_id}: {error_a}")
                continue
            
            is_valid_b, error_b = validate_starting_pdb(structure, chain_b_id, pdb_id)
            if not is_valid_b:
                print(f"Skipping {pdb_id}: {error_b}")
                continue
            
            # Extract sequences
            seq_a = chain2seq(structure, chain_a_id)
            seq_b = chain2seq(structure, chain_b_id)

            # Use helper to isolate each requested chain into its own structure.
            # This avoids warnings and guarantees removing any ligands or extra
            # chains automatically.
            struct_a = copy_structure_with_only_chain(structure, chain_a_id)
            # model id may not be 0 (copy helper uses 1), so grab it dynamically
            first_model = next(struct_a.get_models())
            change_chain_id(struct_a, first_model.id, chain_a_id, 'A')
            
            # Remove HETATM records from template_A before validation
            remove_hetatm_from_structure(struct_a)
            
            # Validate template_A
            is_valid_template_a, error_template_a = validate_template(struct_a, 'A', 'template_A')
            if not is_valid_template_a:
                print(f"Skipping {pdb_id}: {error_template_a}")
                continue
            
            path_a_cif = os.path.join(work_dir, f"template_A.cif")
            save_cif_with_header(struct_a, 'A', path_a_cif)

            struct_b = copy_structure_with_only_chain(structure, chain_b_id)
            first_model = next(struct_b.get_models())
            change_chain_id(struct_b, first_model.id, chain_b_id, 'B')
            
            # Remove HETATM records from template_B before validation
            remove_hetatm_from_structure(struct_b)
            
            # Validate template_B
            is_valid_template_b, error_template_b = validate_template(struct_b, 'B', 'template_B')
            if not is_valid_template_b:
                print(f"Skipping {pdb_id}: {error_template_b}")
                continue
            
            path_b_cif = os.path.join(work_dir, f"template_B.cif")
            save_cif_with_header(struct_b, 'B', path_b_cif)

            # Generate JSONs
            # 1. Binary (Protein A + Protein B)
            binary_json = create_af3_json(
                f"{pdb_id}_binary", seq_a, path_a_cif, seq_b, path_b_cif)
            binary_json_path = os.path.join(work_dir, f"{pdb_id}_binary.json")
            with open(binary_json_path, "w") as jf:
                json.dump(binary_json, jf, indent=2)
            
            # Copy binary JSON to json_bin folder
            binary_json_dest = os.path.join(json_bin_dir, f"{pdb_id}_binary.json")
            shutil.copy(binary_json_path, binary_json_dest)

            # 2. Ternary (Protein A + Protein B + Ligand)
            ternary_json = create_af3_json(
                f"{pdb_id}_ternary", seq_a, path_a_cif, seq_b, path_b_cif, smiles, ligand_id)
            ternary_json_path = os.path.join(work_dir, f"{pdb_id}_ternary.json")
            with open(ternary_json_path, "w") as jf:
                json.dump(ternary_json, jf, indent=2)
            
            # Move ternary JSON to json_ter folder
            ternary_json_dest = os.path.join(json_ter_dir, f"{pdb_id}_ternary.json")
            shutil.copy(ternary_json_path, ternary_json_dest)

            print(f"Processed {pdb_id}: Binary JSON copied to {json_bin_dir}, Ternary JSON copied to {json_ter_dir}")

if __name__ == "__main__":
    # Ensure you have your .cif files in PDB_STORAGE_DIR or use Bio.PDB.PDBList to download them.
    process_complexes(CSV_INPUT)



"""
what to do next:
sbatch /work/lpdi/users/eline/CID_scoring/scripts/run_alphafold.sh -i /work/lpdi/users/eline/CID_scoring/benchmark/af3_inputs/json_bin -o /work/lpdi/users/eline/CID_scoring/benchmark/binary_out --no-msa --num_recycles 3
sbatch /work/lpdi/users/eline/CID_scoring/scripts/run_alphafold.sh -i /work/lpdi/users/eline/CID_scoring/benchmark/af3_inputs/json_ter -o /work/lpdi/users/eline/CID_scoring/benchmark/ternary_out --no-msa --num_recycles 3

sbatch "$SDIR/scripts/run_alphafold.sh" -i "$WDIR/output/af3ternary/json/" -o "$WDIR/output/af3ternary/" --no-msa --num_recycles 3

"""