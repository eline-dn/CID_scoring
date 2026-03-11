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
python "$SDIR/scripts/AF3_output_module.py" --ref_pdb ./input/binder_refs/*.pdb --AF3_outs ./output/pass_af3/pass_pyr/af3_bin/* --outdir ./output/pass_af3/pass_pyr/af3_bin --prefix ""

# pyrosetta scoring
#CONDAPATH="/work/lpdi/users/eline/miniconda3"
python "$SDIR/scripts/2.1_pyrosetta_module.py" --pdb ./output/pass_af3/pass_pyr/af3_bin/*/*_model.pdb --outdir ./output/pass_af3/pass_pyr/af3_bin --prefix ""
echo "done pyrosetta af3 binary"
# ipSAE and dockq
#ids=$(printf '%s\n' ./output/pass_af3/pass_pyr/af3_bin/*/ | sed 's:/$::; s:.*/::')

#python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids --out-csv ./output/pass_af3/pass_pyr/af3_bin/ipsae_and_ipae.csv --af3-dir ./output/pass_af3/pass_pyr/af3_bin/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 
#echo "done ipsae af3 binary"
# dockQ bis:
conda deactivate
python "$SDIR/scripts/pDockQ.py" --native-pdbs input/binder_refs/*.pdb --model-pdbs ./output/pass_af3/pass_pyr/af3_bin/*/*_model.pdb --out-csv ./output/pass_af3/pass_pyr/af3_bin/pdockQ2bin.csv 
echo "done dockq af3 binary"
# plip:
#mkdir output/af3binary/plip
mkdir "output/pass_af3/pass_pyr/af3_bin/plip"
#plip -f output/pass_af3/pass_pyr/af3_bin/t3*/*_model.pdb --inter B  -tx --out "output/pass_af3/pass_pyr/af3_bin/plip"
#python "$SDIR/scripts/3_plip_interaction_profile.py" --xml_files output/pass_af3/pass_pyr/af3_bin/plip/*.xml --outdir output/pass_af3/pass_pyr/af3_bin/plip
echo "done af3 output in  $(($(date +%s) - START_TIME)) seconds"

#/work/lpdi/users/eline/CID/b2_FUN_1Z9Y/input/binder_refs/t3_2_22_1_T0.2_lTp0.2_dcut500.0_seq1.pdb

python "$SDIR/scripts/pDockQ.py" --native-pdbs "/work/lpdi/users/eline/CID/b2_FUN_1Z9Y/input/binder_refs/t3_2_22_1_T0.2_lTp0.2_dcut500.0_seq1.pdb" --model-pdbs "output/missing_binder/t3_2_22_1_t0.2_ltp0.2_dcut500.0_seq1/*_model.pdb" --out-csv "./output/missing_binder/pdockQ2bin.csv"