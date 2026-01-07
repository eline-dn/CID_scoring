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
#WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN_p1/" # change project dir here
WDIR= "$1"
cd "$WDIR"

mkdir output/colabdesign

#-------------- colab design reprediction:----------------------
#BC_env # bc env + module load cuda, etc...
source /work/lpdi/users/mpacesa/Pipelines/miniforge3/bin/activate /work/lpdi/users/mpacesa/Pipelines/miniforge3/envs/BindCraft_kuma ; module load gcc/13.2 ; module load cuda/12.4.1 ; module load cudnn/8.9.7.29-12
python "$SDIR/scripts/Colab_Design_module.py" \
  --ref_pdb input/binder_refs/*.pdb \
  --outdir output/colabdesign \
  --prefix 1Z9Y_FUN

# output pyrosetta scoring:
python "$SDIR/scripts/pyrosetta_module.py" --pdb output/colabdesign/*_model*.pdb --outdir output/colabdesign --prefix colab_pyr
