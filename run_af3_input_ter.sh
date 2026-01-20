#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 20:00:00
#SBATCH --account lpdi
#SBATCH --job-name af3input
#SBATCH --output=logs/af3input_ter_%j.out

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring"

WDIR="$1"
target_id="$2"
target_template="$3"
smiles="$4"
lig_name="$5"

cd "$WDIR"
mkdir "$WDIR/output/af3binary/"
mkdir "$WDIR/output/af3binary/json"
echo "Working directory: $WDIR"

#------------ af3 ternary reprediction:-----------------
# input preparation:
mkdir "$WDIR/output/af3ternary/"
mkdir "$WDIR/output/af3ternary/json"
python "$SDIR/scripts/AF3_input_module.py"  --structure ./input/binder_refs/*.pdb  --target_id "$target_id" --target_template "$target_template"  --smiles "$smiles" --lig_name "$lig_name" --outdir "$WDIR/output/af3ternary/json"
# "c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl"
