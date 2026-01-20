#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 70:00:00
#SBATCH --account lpdi
#SBATCH --job-name af3output_ter
#SBATCH --output=logs/af3output_ter_%j.out

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring"
START_TIME=$(date +%s)
# WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN_p2" # change project dir here
WDIR="$1"
params="$2" # path to the ligand param file
lig_name="$3"

cd "$WDIR"
source ~/.bashrc
#conda activate rosetta_scoring
#------------ af3 ternary reprediction:-----------------
#output analysis
python "$SDIR/scripts/AF3_output_module.py" --ref_pdb ./input/binder_refs/*.pdb --AF3_outs ./output/af3ternary/* --lig_name "$lig_name" --outdir ./output/af3ternary --prefix af3ter

#conda deactivate
# ipSAE and dockq
ids=$(printf '%s\n' ./output/af3ternary/*/ | sed 's:/$::; s:.*/::')

python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids --out-csv ./output/af3ternary/ipsae_and_ipae.csv --af3-dir ./output/af3ternary/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 
echo "done ipsae af3 ternary"
# dockQ bis:
python "$SDIR/scripts/pDockQ.py" --native-pdbs input/binder_refs/*.pdb --model-pdbs ./output/af3ternary/*/*_model.pdb --out-csv ./output/af3ternary/pdockQ2ter.csv --ternary
echo "done pDockQ af3 ternary"

# merging csv and filtering out bad models
python "$SDIR/scripts/merge_af3_ter_csvs.py" --confidence_csv ./output/af3ternary/*_AF3_*reprediction_metrics.csv --pdockq_csv ./output/af3ternary/*pdockQ2ter.csv --ipsae_csv ./output/af3ternary/ipsae_and_ipae.csv --out_dir ./output/af3ternary --filters "./input/filters.json"

# plip:
#mkdir output/af3ternary/plip

#plip -f output/af3ternary/*/*_model.pdb -tx --out output/af3ternary/plip
#python "$SDIR/scripts/plip_interaction_profile.py" --xml_files ./output/af3ternary/plip/*.xml --outdir ./output/af3ternary --prefix af3ter --ternary

echo "done af3 output in  $(($(date +%s) - START_TIME)) seconds"
