# CID_scoring
A set of metrics to score CID and MG complexes


# Usage of each module:

## Reprediction of binary complex with Colab Design:
```
python ./scripts/Colab_Design_module.py \
--ref_pdb structure1_binary_complex.pdb structure2_binary_complex.pdb \
--outdir ./directory \
--prefix prefix \
```

ref_pdb: List of ref PDBs, used to extract sequences, target template, and as references for RMSDs. Ensure the chain id's are correct.
outdir (optional, default is WD) : output folder for csv file and repredictions
prefix (optional, default is WD) : prefix for csv file

## Running AF3 on binary or ternary complexes:
1- Prepare the input json file
```
python ./scripts/AF3_input_module.py \
--pdb structure1_ternary_complex.pdb
--target_id 1Z9Y
--smiles c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl
--lig_name FUN
--outdir "./out"
```
pdb: List of Input PDBs, used to extract sequences
smiles: **(for ternary complexes only)** smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand
lig_name: **(for ternary complexes only)** name used for ligand parametrization, eg  FUN ligand
outdir: output folder for json inputs, optional

2- run AF3
```
sbatch ./scripts/run_alphafold.sh -i input_folder -o output_folder --no-msa --recycles 5
```
input_folder: with json files

3- Analyze the output structure (RMSD, confidences)
```
python ./scripts/AF3_output_module.py \
--ref_pdb ref_structure1_ternary_complex.pdb \
--AF3_outs output_folder/* \
--lig_name FUN \
--outdir . \
--prefix pre\
```
ref_pdb: List of ref PDBs, used to extract sequences, as target template, and as references for RMSDs. Ensure the id is the same as the folder with repredicted strucures and confidences. Also ensure the chain id's are matching.
AF3_outs: List of AF3 output folders. Ensure the id is the same as the ref_pdb.
lig_name: name used for ligand parametrization, eg  FUN ligand
outdir: output folder for csv file"
prefix: prefix for csv file


## PyRosetta Scoring (inspired from BindCraft)
```
python ./scripts/pyrosetta_module.py \
--pdb repredicted_structure1_ternary_complex.pdb \
--mk_params \
--smiles c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl \
--lig_name FUN \
--outdir . \
--prefix pre\
```
pdb: List of Input PDBs
mk_params: need to create a Params files in the pdb files'dir? Set to true if the pdbs are ternary complexes. in that case, add a smiles for the ligand
smiles: **(for ternary complexes only)** smiles used for ligand parametrization, eg c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl for the FUN ligand
lig_name: **(for ternary complexes only)** name used for ligand parametrization, eg  FUN ligand
outdir: output folder for csv file"
prefix: prefix for csv file


## Compute DockQ (https://github.com/DigBioLab/de_novo_binder_scoring/tree/main)
```
python ./scripts/dockQ.py \
  --run-csv ./example_output/run.csv \
  --input-pdbs ./example_output/input_pdbs/ \
  --folder af3:./example_output/AF3/pdbs/ \
  --out-csv ./example_output/dockQ.csv 
```

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
