# General workflow:

Prepare data_complexes.csv

```
AF3_input_prep


sbatch /work/lpdi/users/eline/CID_scoring/scripts/run_alphafold.sh -i /work/lpdi/users/eline/CID_scoring/benchmark/af3_inputs/json_bin -o /work/lpdi/users/eline/CID_scoring/benchmark/binary_out --no-msa --num_recycles 3
sbatch /work/lpdi/users/eline/CID_scoring/scripts/run_alphafold.sh -i /work/lpdi/users/eline/CID_scoring/benchmark/af3_inputs/json_ter -o /work/lpdi/users/eline/CID_scoring/benchmark/ternary_out --no-msa --num_recycles 3
```
# next: AF3 output analysis:
We run classical confidence scores 

```bash
#confidence ter and bin
# need to modify script AF3_output_module.py




# pyrosetta scoring bin
#CONDAPATH="/work/lpdi/users/eline/miniconda3"
python "$SDIR/scripts/2.1_pyrosetta_module.py" --pdb ./binary_out/*/*_model.cif --outdir ./ --prefix "bin_" --cif


``` 
and ipsae script + pdockQ

# ipSAE and dockq
```
SDIR="/work/lpdi/users/eline/CID_scoring"
ids_ter=$(printf '%s\n' ./ternary_out/*/ | sed 's:/$::; s:.*/::')
ids_bin=$(printf '%s\n' ./binary_out/*/ | sed 's:/$::; s:.*/::')

python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids_ter --out-csv ./ter_ipsae_and_ipae.csv --af3-dir ./ternary_out/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 

python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids_bin --out-csv ./bin_ipsae_and_ipae.csv --af3-dir ./binary_out/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 

# dockQ bis:
python "$SDIR/scripts/pDockQ.py" --native-pdbs pdb_files/*.cif --model-pdbs ./ternary_out/*/*_model.pdb --out-csv ./pdockQ2ter.csv --ternary

python "$SDIR/scripts/pDockQ.py" --native-pdbs pdb_files/*.cif --model-pdbs ./binary_out/*/*_model.pdb --out-csv ./pdockQ2ter.csv --ternary
```