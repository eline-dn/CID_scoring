"""
goal: input: a list of complexes, binary or ternary
(repredict complexes
extract confidence metrics
compute rmsds) in part 2
"""

# load dependencies
import os,  sys
import numpy as np
import argparse
import pandas as pd
import time

# import AF3_utils.py
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/../functions")
from AF3_utils import *
from biopython_utils import *


# general template json file:
template_json_path=f"{SCRIPT_PATH}/no_msa_template.json"
script_start_time = time.time()







######## 

parser = argparse.ArgumentParser()
parser.add_argument("--structure", nargs="+", type=str, help="List of Input CIF or pdb, used to extract sequences and target template. Target must be chain A and the first residues in the structure") # list of pdb files 
parser.add_argument("--target_id", required=True, type=str, help="target id")
parser.add_argument("--target_template", required=True, type=str, help="target template to input to af3, must be in a CIF format with a release date before the train/test cutoff of AF3. Give an absolute path")
parser.add_argument("--smiles", required=False, type=str, help="smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand")
parser.add_argument("--lig_name", required=False, type=str, help="name used for ligand parametrization, eg  FUN ligand")
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file")
#parser.add_argument("--prefix", required=False, type=str, help="prefix for csv file")
args = parser.parse_args()

if ((args.smiles is not None) and (args.lig_name is None)) or ((args.smiles is None) and (args.lig_name is not None)):
  raise ValueError(" if you go for a ternary complex reprediction, please specify smiles AND ligand name")

if (args.smiles is not None) and (args.lig_name is not None):
  smiles=args.smiles
  lig_name=args.lig_name
  ternary=True
  print("reprediction of ternary complexes, ensure the given pdb look like: A target, B binder, ligand")
else:
  ternary=False
  
outdir=args.outdir

""" to write a custom target template structure:
# template must be in a cif format for AF3, converting if given as a pdb
# also: template structure must be filtered to a single chain.
if '.cif' in args.structure[0]:
  format="CIF"
  struct=load_CIF(args.structure[0])
  struct=copy_structure_with_only_chain(struct, "A")
  if not args.structure[0].startswith("/"): # if an absolute path is not provided:
    # template absolute path:
    target_template=os.getcwd()
    target_template+=f"/{args.structure[0]}"
    target_template=target_template.replace("//", "/")
  else:
    target_template=args.structure[0]
  write_cif(struct,target_template)
  
elif '.pdb' in args.structure[0]:
  format="PDB"
  struct=load_PDB(args.structure[0])
  struct=copy_structure_with_only_chain(struct, "A")
  # path to save cif tempalte:
  if not outdir.startswith("/"): # if outdir is not an absolute path
    target_template=os.getcwd()
    target_template+=f"/{outdir}/template_{args.target_id}.cif"
    target_template=target_template.replace("//", "/")
  else:
    target_template=f"{outdir}/template_{args.target_id}.cif"
  write_cif(struct,target_template)
  
else: 
  raise ValueError("input structure format must be cif or pdb")
"""
target_template=args.target_template
# prepare target specific json input file
if ternary:
  target_json=prep_no_msa_target_input_json(gen_template_path=template_json_path, outdir=outdir, target_id=args.target_id, target_cif=target_template, target_chain_in_cif="A", smiles=smiles, ligand_id=lig_name)
else:
  target_json=prep_no_msa_target_input_json(gen_template_path=template_json_path, outdir=outdir, target_id=args.target_id, target_cif=target_template, target_chain_in_cif="A", smiles=None, ligand_id=None)
  
##### prepare inputs in the same folder, everything repredicted in one aF3 command
for struct in args.structure:
  if format=='CIF':
    id=os.path.basename(struct).replace(".cif", "")
    structure=load_CIF(struct)
  else:
    id=os.path.basename(struct).replace(".pdb", "")
    structure=load_PDB(struct)
  seq_binder=chain2seq(structure, "B", make_str=True)
  json=prep_binder_input_json(template_path=target_json, outdir=outdir ,id=id ,seq_binder=seq_binder)


elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
n_binder=len(args.structure)
print(f"Finished AF3 input preparation for {n_binder} complexes. Script execution took: "+elapsed_text)

"""
what to do next:
sbatch ./run_alphafold.sh -i input2/1 -o output2 --no-msa
sbatch ./run_alphafold.sh -i input2/2 -o output2 --no-msa
sbatch ./run_alphafold.sh -i input2/3 -o output2 --no-msa


sbatch run_alphafold.sh -i /work/lpdi/users/dobbelst/tools/alphafold3_examples/af_input/fold_input_singleseq.json -o <OUTPUT_DIR> --no-msa
"""
