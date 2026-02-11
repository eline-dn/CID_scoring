#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 72:00:00
#SBATCH --account lpdi
#SBATCH --job-name plip
#SBATCH --output=logs/plip_%j.out

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring"
START_TIME=$(date +%s)
# WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN_p2" # change project dir here
WDIR="$1"

cd "$WDIR"
source ~/.bashrc

# plip:
mkdir output/pass_af3/pass_pyr/plip

plip -f output/pass_af3/pass_pyr/*/*_model.pdb -tx --out output/pass_af3/pass_pyr/plip
python "$SDIR/scripts/3_plip_interaction_profile.py" --xml_files ./output/pass_af3/pass_pyr/plip/*.xml --outdir ./output/pass_af3/pass_pyr --ternary

if [ ! -f "./output/pass_af3/pass_pyr/BL__plip_interactions.csv" ] ; then
    echo "plip profiling step failed"
    exit 1
fi

# merging the plip scoring with the previous scoring metrics and filtering again

python "$SDIR/scripts/3.2_merge_plip_ter.py" --plip_csv ./output/pass_af3/pass_pyr/BL__plip_interactions.csv --accepted1_csv ./output/pass_af3/pass_pyr/pass_pyr_ter.csv --out_dir ./output/pass_af3/pass_pyr

echo "Number of binders remaining in hbond folder:"
ls ./output/pass_af3/pass_pyr/hbond/ | wc -l

if [ ! -f ".output/pass_af3/pass_pyr/all_af3_pyr_plip_ter.csv" ] ; then
    echo "plip filtering step failed"
    exit 1
fi

echo "done plip scoring in  $(($(date +%s) - START_TIME)) seconds"
