"""
usage: 
python "$SDIR/scripts/merge_af3_ter_csvs.py" --confidence_csv ./output/af3ternary/af3_ter_pyr_scoring.csv --pdockq_csv ./output/af3ternary/pdockQ2ter.csv --ipsae_csv ./output/af3ternary/ipsae_and_ipae.csv --out_csv ./output/af3ternary/merged_af3_ter_scoring.csv --filters "./input/filters.json"
Collect the informations in ./output/af3ternary/accepted2/accepted_metrics2.csv and in .output/af3binary/*.csv (af3 confidences, ipae/ipSAE, pdockQ, pyrosetta)

"""

import pandas as pd
import argparse
import json 
# load dependencies
import os,  sys
import numpy as np
import time
import glob
import re

# --- HELPER FUNCTIONS ---
def clean_id(id_str):
    """
    Clean structure ID by removing the '_model' suffix and trailing digits.
    Example:
        t2_2_7_3_t0.1_ltp0.2_dcut15.0_model0 -> t2_2_7_3_t0.1_ltp0.2_dcut15.0
    """
    return re.sub(r'_model\d*$', '', id_str).lower()


def clean(df: pd.DataFrame, folder_prefix) -> pd.DataFrame:
    """
    Collapse rows whose ID ends with '_modelX' (X = single digit)
    by averaging numeric columns.
    clean id
    add prefix
    """
    df = df.copy()
    # Identify the ID column: common possibilities
    id_col_candidates = [c for c in df.columns if c.lower() in ["id", "binder_id","ter_id", "ID"]]
    if not id_col_candidates:
        raise ValueError("No ID column found in dataframe")
    id_col = id_col_candidates[0]
    # Standardize the ID column
    df = df.rename(columns={id_col: "ID"})
  
    # clean ID
    df["ID"] = df["ID"].astype(str)
    df["ID_clean"] = df["ID"].apply(clean_id)

    # identify numeric columns only
    numeric_cols = df.select_dtypes(include="number").columns.tolist()

    # group by cleaned ID and average numeric columns
    grouped = (
        df.groupby("ID_clean", as_index=False)[numeric_cols]
          .mean()
    )

    # rename back to ID
    grouped = grouped.rename(columns={"ID_clean": "ID"})
    # Prefix other columns
    other_cols = [c for c in grouped.columns if c != "ID"]
    grouped = grouped.rename(columns={c: f"{folder_prefix}_{c}" for c in other_cols})
    return grouped


#----------------parse arguments----------------
parser = argparse.ArgumentParser()
parser.add_argument("--confidence_csv", required=True, type=str, help="csv file with af3 confidence metrics")
parser.add_argument("--pdockq_csv", required=True, type=str, help="csv file with pdockQ2 scoring metrics")
parser.add_argument("--ipsae_csv", required=True, type=str, help="csv file with ipsae and ipae scoring metrics")
parser.add_argument("--pyr_csv", required=True, type=str, help="csv file with ipsae and ipae scoring metrics")
parser.add_argument("--out_dir", required=True, type=str, help="output file for merged csv and good/bad models subfolders")
parser.add_argument("--ter_merged", required=True, type=str, help="csv file with af3 ternary information for our selected binders")
args = parser.parse_args()

# merge the csvs on the id column after loading and cleaning them
script_start_time = time.time()
out_dir = args.out_dir

conf_df = pd.read_csv(args.confidence_csv)
pdockq_df = pd.read_csv(args.pdockq_csv)
ipsae_df = pd.read_csv(args.ipsae_csv)
pyr_df = pd.read_csv(args.pyr_csv)
ter_merged_df = pd.read_csv(args.ter_merged)


# ---------- Merge ---
folder_dfs = []
merged_df = None
ter_merged_df = clean(ter_merged_df, folder_prefix="ter_")
folder_dfs.append(ter_merged_df)
for df in [conf_df, pdockq_df, ipsae_df, pyr_df]:
    #df = clean_dataframe(df, folder_prefix)
    df = clean(df, folder_prefix="bin_")
    folder_dfs.append(df)
# Merge CSVs on "ID"
merged_df = folder_dfs[0]
for df in folder_dfs[1:]:
    merged_df = merged_df.merge(df, on="ID", how="inner")
# rule: metrics from the ternary complexes start with "ter_" and metrics from the binary complexes start with "bin_
"
# ----------- Processing the "diff" columns-------
"""diff_pairs = {
    "binder_score":(),
    "surface_hydrophobicity": () interface_sc interface_dG interface_dSASA _interface_dG_SASA_ratio 
    _interface_nres _interface_hbond_percentage

    ter_iptm

    ter_d_binding_site  _ipSAE_min  af3_LIS  
}"""
# col ends with: 
metrics_to_diff = ["binder_score", "surface_hydrophobicity", "interface_sc",
 "interface_dG", "interface_dSASA", "interface_dG_SASA_ratio",
    "_interface_nres", "_interface_hbond_percentage",
    "ter_iptm","ter_d_binding_site",  "_ipSAE_min", "af3_LIS"]
"""
def diff_metric(row, metric):
    ter_cols = [col for col in merged_df.columns if col.startswith("ter_") and metric in col]
    bin_cols= [col for col in merged_df.columns if col.startswith("bin_") and metric in col]
    print(f"Metric: {metric}, ter_cols: {ter_cols}, bin_cols: {bin_cols}")
    if len(ter_cols)> 1 or len(bin_cols)> 1:
        print(f"Warning: multiple columns found for metric {metric}, taking the first one")
    if len(ter_cols)< 1 or len(bin_cols)< 1:
        raise ValueError(f"Warning: no column found for metric {metric}")
    diff= row[str(ter_cols[0])]- row[str(bin_cols[0])]
    return diff

for metric in metrics_to_diff:
    merged_df[f"delta_{metric}"]=merged_df.apply(diff_metric, metric=metric)
"""
# --- SAVE ---
if merged_df is not None:
output_file = os.path.join("/work/lpdi/users/eline/CID/FUN_1Z97/output", "merged_full.csv")
merged_df.to_csv(output_file, index=False)
    print(f"Metrics saved to {output_file}")
else:
    print("No CSV files found.")



