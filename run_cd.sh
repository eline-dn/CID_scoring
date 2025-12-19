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
SDIR="/work/lpdi/users/eline/CID_scoring/"
WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN_p2/"

cd "$WDIR"

mkdir output/colabdesign

#-------------- colab design reprediction:----------------------
BC_env
python "$SDIR/scripts/Colab_Design_module.py" \
  --ref_pdb input/binder_refs/*.pdb \
  --outdir output/colabdesign \
  --prefix 1Z9Y_FUN

# output pyrosetta scoring:
python "$SDIR/scripts/pyrosetta_module.py" --pdb output/colabdesign/*_model*.pdb --outdir output/colabdesign --prefix colab_pyr
