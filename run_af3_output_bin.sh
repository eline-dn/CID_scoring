#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 70:00:00
#SBATCH --account lpdi
#SBATCH --job-name af3output_bin
#SBATCH --output=logs/af3output_bin_%j.out

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring"
START_TIME=$(date +%s)
WDIR="$1"

cd "$WDIR"
source ~/.bashrc

#------------ af3 binary reprediction:-----------------
#output analysis
conda activate rosetta_scoring
#python "$SDIR/scripts/AF3_output_module.py" --ref_pdb ./input/binder_refs/*.pdb --AF3_outs ./output/af3binary/* --outdir ./output/af3binary --prefix af3bin

# pyrosetta scoring
#CONDAPATH="/work/lpdi/users/eline/miniconda3"
python "$SDIR/scripts/pyrosetta_module.py" --pdb ./output/af3binary/*/*_model.pdb --outdir ./output/af3binary --prefix af3_bin_pyr
echo "done pyrosetta af3 binary"
# ipSAE and dockq
#ids=$(printf '%s\n' ./output/af3binary/*/ | sed 's:/$::; s:.*/::')

#python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids --out-csv ./output/af3binary/ipsae_and_ipae.csv --af3-dir ./output/af3binary/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 
echo "done ipsae af3 binary"
# dockQ bis:
conda deactivate
#python "$SDIR/scripts/pDockQ.py" --native-pdbs input/binder_refs/*.pdb --model-pdbs ./output/af3binary/*/*_model.pdb --out-csv ./output/af3binary/pdockQ2bin.csv 
echo "done dockq af3 binary"
# plip:
#mkdir output/af3binary/plip

#plip -f "$WDIR/output/af3binary/*/*_model.pdb" --inter B  -tx --out "$WDIR/output/af3binary/plip"
#python "$SDIR/scripts/plip_interaction_profile.py" --xml_files output/af3binary/plip/*.xml --outdir output/af3binary/plip --prefix af3bin
echo "done af3 output in  $(($(date +%s) - START_TIME)) seconds"
