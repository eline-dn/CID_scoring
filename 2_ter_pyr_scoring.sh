#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 72:00:00
#SBATCH --account lpdi
#SBATCH --job-name af3output_ter
#SBATCH --output=logs/af3output_ter_bis_%j.out

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring"
START_TIME=$(date +%s)
# WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN_p2" # change project dir here
WDIR="$1"
params="$2" # path to the ligand param file
lig_name="$3"

cd "$WDIR"
source ~/.bashrc

# plip:
#mkdir output/af3ternary/plip

#plip -f output/af3ternary/accepted3/*/*_model.pdb -tx --out output/af3ternary/plip
#python "$SDIR/scripts/plip_interaction_profile.py" --xml_files ./output/af3ternary/plip/*.xml --outdir ./output/af3ternary --prefix af3ter --ternary

# running the pyrosetta relaxation and scoring on the succesful binders:
# pyrosetta scoring
conda activate rosetta_scoring
CONDAPATH="/work/lpdi/users/eline/miniconda3"
"$CONDAPATH/envs/rosetta_scoring/bin/python" "$SDIR/scripts/2.1_pyrosetta_module.py" --pdb ./output/pass_af3/*/*_model.pdb --params "$params" --lig_name "$lig_name" --outdir ./output/pass_af3 --prefix ""

if [! -f ".output/pass_af3/_ter_pyrosetta_metrics.csv"]; then 
    echo "rosetta scoring step failed"
    exit 1
fi
echo "done pyrosetta af3 ternary"
#conda deactivate

# merging the pyrosetta scoring with the previous scoring metrics and filtering again
python "$SDIR/scripts/2.2_merge_pyr_ter.py" --pyr_csv ./output/pass_af3/*pyrosetta_metrics.csv --accepted1_csv ./output/pass_af3/pass_af3_ter.csv --out_dir ./output/pass_af3 --filters "./input/pyr_filters.json"

echo "Number of binders remaining in pass_pyr folder:"
ls ./output/pass_af3/pass_pyr/ | wc -l

if [! -f ".output/pass_af3/all_af3_pyr_ter.csv"]; then
    echo "filtering step failed"
    exit 1
fi

echo "done pyr scoring in  $(($(date +%s) - START_TIME)) seconds"
