# setup working directory
SDIR="/work/lpdi/users/eline/CID_scoring/"
WDIR="/work/lpdi/users/eline/CID/1Z9Y_FUN"

cd "$WDIR"
# conda env:
CONDAPATH="/work/lpdi/users/eline/miniconda3"  # edit this depending on where your Conda environments live
# "$CONDAPATH/envs/rosetta_scoring/bin/python"
/work/lpdi/users/eline/miniconda3/envs/rosetta_scoring/bin/python

mkdir input
mkdir output
mkdir input/binder_refs
mkdir output/colabdesign
mkdir output/af3binary
mkdir output/af3binary/json
mkdir output/af3ternary
mkdir output/af3ternary/json

#-------------- colab design reprediction:----------------------
BC_env
python "$SDIR/scripts/Colab_Design_module.py" \
  --ref_pdb input/binder_refs/t2_2_7_*.pdb \
  --outdir output/colabdesign \
  --prefix 1Z9Y_FUN
# output pyrosetta scoring:
python "$SDIR/scripts/pyrosetta_module.py" --pdb output/colabdesign/*_model*.pdb --outdir output/colabdesign --prefix colab_pyr


#------------ af3 ternary reprediction:-----------------
# input preparation:
python "$SDIR/scripts/AF3_input_module.py"  --structure ./input/binder_refs/t2_2_7_*.pdb  --target_id 1Z9Y --target_template "$WDIR/input/1Z9Y_target_template.cif"  --smiles "c1cc(oc1)CNc2cc(c(cc2C(=O)O)S(=O)(=O)N)Cl" --lig_name FUN --outdir "$WDIR/output/af3ternary/json"

# run af3
sbatch "$SDIR/scripts/run_alphafold.sh" -i "$WDIR/output/af3ternary/json/" -o "$WDIR/output/af3ternary/" --no-msa --num_recycles 5

#output analysis
python "$SDIR/scripts/AF3_output_module.py" --ref_pdb ./input/binder_refs/t2_2_7_*.pdb --AF3_outs ./output/af3ternary/* --lig_name FUN --outdir ./output/af3ternary --prefix af3ter

# pyrosetta scoring
#CONDAPATH="/work/lpdi/users/eline/miniconda3"
#"$CONDAPATH/envs/rosetta_scoring/bin/python" "$SDIR/scripts/pyrosetta_module.py" --pdb ./output/af3ternary/*/*_model.pdb --params "$WDIR/input/FUN.params" --lig_name FUN --outdir ./output/af3ternary --prefix af3_ter_pyr

# ipSAE and dockq
ids=$(printf '%s\n' ./output/af3ternary/t2*/ | sed 's:/$::; s:.*/::')

python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids --out-csv ./output/af3ternary/ipsae_and_ipae.csv --af3-dir ./output/af3ternary/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 

# dockQ bis:
python "$SDIR/scripts/pDockQ.py" --native-pdbs input/binder_refs/t2_2_7_*.pdb --model-pdbs ./output/af3ternary/*/*_model.pdb --out-csv ./output/af3ternary/pdockQ2ter.csv --ternary

# plip:
mkdir output/af3ternary/plip

plip -f output/af3ternary/*/*_model.pdb -tx --out output/af3ternary/plip
python "$SDIR/scripts/plip_interaction_profile.py" --xml_files ./output/af3ternary/plip/*.xml --outdir ./output/af3ternary --prefix af3ter --ternary


#------------ af3 binary reprediction:-----------------
# input preparation:
python "$SDIR/scripts/AF3_input_module.py"  --structure ./input/binder_refs/t2_2_7_*.pdb  --target_id 1Z9Y --target_template "$WDIR/input/1Z9Y_target_template.cif"  --outdir "$WDIR/output/af3binary/json"

# run af3
sbatch "$SDIR/scripts/run_alphafold.sh" -i "$WDIR/output/af3binary/json/" -o "$WDIR/output/af3binary/" --no-msa --num_recycles 5

#output analysis
python "$SDIR/scripts/AF3_output_module.py" --ref_pdb ./input/binder_refs/t2_2_7_*.pdb --AF3_outs ./output/af3binary/* --outdir ./output/af3binary --prefix af3bin

# pyrosetta scoring
CONDAPATH="/work/lpdi/users/eline/miniconda3"
"$CONDAPATH/envs/rosetta_scoring/bin/python" "$SDIR/scripts/pyrosetta_module.py" --pdb ./output/af3binary/*/*_model.pdb --outdir ./output/af3binary --prefix af3_bin_pyr

# ipSAE and dockq
ids=$(printf '%s\n' ./output/af3binary/t2*/ | sed 's:/$::; s:.*/::')

python "$SDIR/scripts/run_ipsae_batch.py" --id_list $ids --out-csv ./output/af3binary/ipsae_and_ipae.csv --af3-dir ./output/af3binary/ --ipsae-script-path "$SDIR/functions/ipsae_w_ipae.py" --specific-chainpair-ipsae "A:B,B:A"  --pae-cutoff 10 --overwrite-ipsae --dist-cutoff 10 

# dockQ bis:
python "$SDIR/scripts/pDockQ.py" --native-pdbs input/binder_refs/t2_2_7_*.pdb --model-pdbs ./output/af3binary/*/*_model.pdb --out-csv ./output/af3binary/pdockQ2bin.csv 

# plip:
mkdir output/af3binary/plip

plip -f output/af3binary/*/*_model.pdb --inter B  -tx --out output/af3binary/plip
python "$SDIR/scripts/plip_interaction_profile.py" --xml_files output/af3binary/plip/*.xml --outdir output/af3binary/plip --prefix af3bin
