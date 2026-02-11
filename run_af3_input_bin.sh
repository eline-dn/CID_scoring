#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 20:00:00
#SBATCH --account lpdi
#SBATCH --job-name af3input_bin
#SBATCH --output=logs/af3input_bin_%j.out

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring"

WDIR="$1"
target_id="$2"
target_template="$3"


cd "$WDIR"
#mkdir "$WDIR/output/af3binary_bis/"
mkdir "$WDIR/output/pass_af3/pass_pyr/af3_bin/json"
echo "Working directory: $WDIR"
#------------ af3 binary reprediction:-----------------
# input preparation:
# getting the ids from the accepted2 binders:
shopt -s nocasematch
dirs=$(printf '%s\n' ./output/pass_af3/pass_pyr/hbond/*/ | sed 's:/$::; s:.*/::')
paths=()
for dir in $dirs; do
    while IFS= read -r -d '' file; do
        paths+=("$file")
    done < <(find "$WDIR/input/binder_refs" -maxdepth 1 -type f -iname "${dir}.pdb" -print0)
done


# preparing the inputs:
python "$SDIR/scripts/AF3_input_module.py"  --structure "${paths[@]}"  --target_id "$target_id" --target_template "$target_template"  --outdir "$WDIR/output/pass_af3/pass_pyr/af3_bin/json"

