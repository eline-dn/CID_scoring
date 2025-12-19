#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 45G
#SBATCH --gres gpu:1
#SBATCH --partition h100
#SBATCH --time 10:00:00
#SBATCH --account lpdi
#SBATCH --job-name cd2

# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring/"
WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN_p2/" # change project dir here

cd "$WDIR"

#------------ af3 binary reprediction:-----------------
# input preparation:
python "$SDIR/scripts/AF3_input_module.py"  --structure ./input/binder_refs/*.pdb  --target_id 1Z9Y --target_template "$WDIR/input/1Z9Y_target_template.cif"  --outdir "$WDIR/output/af3binary/json"


#------------ af3 ternary reprediction:-----------------
# input preparation:
python "$SDIR/scripts/AF3_input_module.py"  --structure ./input/binder_refs/*.pdb  --target_id 1Z9Y --target_template "$WDIR/input/1Z9Y_target_template.cif"  --smiles "c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl" --lig_name FUN --outdir "$WDIR/output/af3ternary/json"
