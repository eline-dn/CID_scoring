# CID_scoring
A set of metrics to score CID and MG complexes


# Usage of each module:

## Reprediction of binary complex with Colab Design:
```
python ./scripts/Colab_Design_module.py \
--ref_pdb structure1_binary_complex.pdb structure2_binary_complex.pdb \
--outdir . \
--prefix prefix
```

ref_pdb: List of ref PDBs, used to extract sequences, target template, and as references for RMSDs. Ensure the chain id's are correct.

outdir (optional, default is WD) : output folder for csv file and repredictions

prefix (optional, default is WD) : prefix for csv file

## Running AF3 on binary or ternary complexes:
1- Prepare the input json file
```
 python ./scripts/AF3_input_module.py --structure af3_ter_out.pdb --target_id 1Z9Y --target_template /work/lpdi/users/eline/CID_scoring/af3/input/1Z9Y_clean.cif --smiles "c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl" --lig_name FUN --outdir af3/input
```
pdb: List of Input PDBs, used to extract sequences

target_template: target template to input to af3, must be in a CIF format with a release date before the train/test cutoff of AF3. Give an absolute path

smiles: **(for ternary complexes only)** smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand

lig_name: **(for ternary complexes only)** name used for ligand parametrization, eg  FUN ligand

outdir: output folder for json inputs, optional

tips for the target template file: beware of DOS file format, not parsed by AF3. check eol with `file 1Z9Y_clean.cif` and remove CRLF line terminators with `sed -i 's/\r$//' 1Z9Y_clean.cif`
if your template doesn't have a release date, add
```
#
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
in the cif file's header
```


2- run AF3
```
sbatch ./scripts/run_alphafold.sh -i input_folder -o output_folder --no-msa --recycles 5
```
input_folder: with json files


3- Analyze the output structure (RMSD, confidences)
```
python ./scripts/AF3_output_module.py --ref_pdb af3_ter_out.pdb --AF3_outs /work/lpdi/users/eline/CID_scoring/af3_out/* --lig_name FUN --outdir /work/lpdi/users/eline/CID_scoring/af3_out --prefix pre
```
ref_pdb: List of ref PDBs, used to extract sequences, as target template, and as references for RMSDs. Ensure the id is the same as the folder with repredicted strucures and confidences. Also ensure the chain id's are matching.

AF3_outs: List of AF3 output folders. Ensure the id is the same as the ref_pdb.

lig_name: name used for ligand parametrization, eg  FUN ligand

outdir: output folder for csv file"

prefix: prefix for csv file


## PyRosetta Scoring (inspired from BindCraft)
```
python ./scripts/pyrosetta_module.py --pdb /work/lpdi/users/eline/smol_binder_diffusion_
pipeline/1Z9Yout/4_af3/output/t2_1_100_3_t0.1_ltp0.3_dcut8.0/t2_1_100_3_t0.1_ltp0.3_dcut8.0_model.cif --mk_params --smiles "c1cc(oc1)
CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl" --lig_name FUN --outdir . --prefix pre
```
pdb: List of Input PDBs

mk_params: need to create a Params files in the pdb files'dir? Set to true if the pdbs are ternary complexes. in that case, add a smiles for the ligand

smiles: **(for ternary complexes only)** smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand

lig_name: **(for ternary complexes only)** name used for ligand parametrization, eg  FUN ligand

outdir: output folder for csv file"

prefix: prefix for csv file

NB: run with different csv paths/prefix for binary vs ternary complex bc of column name mismatch!


## Compute DockQ (https://github.com/DigBioLab/de_novo_binder_scoring/tree/main)
```
python ./scripts/dockQ.py --input-pdbs ./refpdbs/ --folder af3:./af3_out/n9_l106_s297124_mpnn1_model1_model0/ --out-csv ./dockQ.csv --verbose
```

input-pdbs: Path to REQUIRED reference/input PDBs (binder_id.pdb).

folder: Repeatable. Format "name:path". Example: --folder af3:/path/to/AF3/pdbs


## Compute ipSAE and interface confidence metrics (https://github.com/DigBioLab/de_novo_binder_scoring/tree/main)
```
python ./scripts/run_ipsae_batch.py \
  --run-csv ./example_output/run.csv \
  --out-csv ./example_output/ipsae_and_ipae.csv \
  --af3-dir ./example_output/AF3 \
  --ipsae-script-path ./functions/ipsae_w_ipae.py

```
There is a possibility to extract specific chain pair ipSAE values: use the argument
```
--specific-chainpair-ipsae "A:D,A:B,A:C" which takes a string formated like the above.
```
It is also possible to specify several thressholds for the AF3 contact prob metrics ``` --confidence-threshold "0.5,0.6,0.7,0.8,0.9" \ ```
