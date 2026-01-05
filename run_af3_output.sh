#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 70:00:00
#SBATCH --account lpdi
#SBATCH --job-name cd2

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring"
WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN_p2" # change project dir here

cd "$WDIR"

#------------ af3 ternary reprediction:-----------------
#output analysis
python "$SDIR/scripts/AF3_output_module.py" --ref_pdb ./input/binder_refs/*.pdb --AF3_outs ./output/af3ternary/* --lig_name FUN --outdir ./output/af3ternary --prefix af3ter

# pyrosetta scoring
CONDAPATH="/work/lpdi/users/eline/miniconda3"
"$CONDAPATH/envs/rosetta_scoring/bin/python" "$SDIR/scripts/pyrosetta_module.py" --pdb ./output/af3ternary/*/*_model.pdb --params "$WDIR/input/FUN.params" --lig_name FUN --outdir ./output/af3ternary --prefix af3_ter_pyr

# ipSAE and dockq
ids=$(printf '%s\n' ./output/af3ternary/t2*/ | sed 's:/$::; s:.*/::')

python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids --out-csv ./output/af3ternary/ipsae_and_ipae.csv --af3-dir ./output/af3ternary/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 

# dockQ bis:
python "$SDIR/scripts/pDockQ.py" --native-pdbs input/binder_refs/*.pdb --model-pdbs ./output/af3ternary/*/*_model.pdb --out-csv ./output/af3ternary/pdockQ2ter.csv --ternary

# plip:
mkdir output/af3ternary/plip

plip -f output/af3ternary/*/*_model.pdb -tx --out output/af3ternary/plip
python "$SDIR/scripts/plip_interaction_profile.py" --xml_files ./output/af3ternary/plip/*.xml --outdir ./output/af3ternary --prefix af3ter --ternary


#------------ af3 binary reprediction:-----------------
#output analysis
python "$SDIR/scripts/AF3_output_module.py" --ref_pdb ./input/binder_refs/*.pdb --AF3_outs ./output/af3binary/* --outdir ./output/af3binary --prefix af3bin

# pyrosetta scoring
CONDAPATH="/work/lpdi/users/eline/miniconda3"
"$CONDAPATH/envs/rosetta_scoring/bin/python" "$SDIR/scripts/pyrosetta_module.py" --pdb ./output/af3binary/*/*_model.pdb --outdir ./output/af3binary --prefix af3_bin_pyr

# ipSAE and dockq
ids=$(printf '%s\n' ./output/af3binary/t2*/ | sed 's:/$::; s:.*/::')

python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids --out-csv ./output/af3binary/ipsae_and_ipae.csv --af3-dir ./output/af3binary/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 

# dockQ bis:
python "$SDIR/scripts/pDockQ.py" --native-pdbs input/binder_refs/*.pdb --model-pdbs ./output/af3binary/*/*_model.pdb --out-csv ./output/af3binary/pdockQ2bin.csv 

# plip:
mkdir output/af3binary/plip

plip -f "$WDIR/output/af3binary/*/*_model.pdb" --inter B  -tx --out "$WDIR/output/af3binary/plip"
python "$SDIR/scripts/plip_interaction_profile.py" --xml_files output/af3binary/plip/*.xml --outdir output/af3binary/plip --prefix af3bin

