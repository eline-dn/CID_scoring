"""
usage: 
python "$SDIR/scripts/merge_af3_ter_csvs.py" --confidence_csv ./output/af3ternary/af3_ter_pyr_scoring.csv --pdockq_csv ./output/af3ternary/pdockQ2ter.csv --ipsae_csv ./output/af3ternary/ipsae_and_ipae.csv --out_csv ./output/af3ternary/merged_af3_ter_scoring.csv --filters "./input/filters.json"
goal: merge the different csv files with scoring metrics for ternary complexes and filter out bad models based on given criteria
put the scoring thresholds in a json file

merge csvs
filter out : plddt > 0.90, iptM, ipSAE, full rmsd and rmsd to binding site, also if no contacts (min 3) <<<<
mv unsuccessful binders (whole af3 folder) to failed folder in af3ter/failed/id2/_model.cif
write metrics to failed1.csv or accepted1.csv

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
parser.add_argument("--out_dir", required=True, type=str, help="output file for merged csv and good/bad models subfolders")
parser.add_argument("--filters", required=True, type=str, help="json file with filtering criteria")
args = parser.parse_args()

# merge the csvs on the id column after loading and cleaning them
script_start_time = time.time()
out_dir = args.out_dir

conf_df = pd.read_csv(args.confidence_csv)
pdockq_df = pd.read_csv(args.pdockq_csv)
ipsae_df = pd.read_csv(args.ipsae_csv)


# ---------- Merge ---
merged_df = None
for df in [conf_df, pdockq_df, ipsae_df]:
    #df = clean_dataframe(df, folder_prefix)
    df = clean(df, folder_prefix="af3_ter")
    folder_dfs.append(df)
# Merge CSVs within the folder on "ID"
if folder_dfs:
    folder_merged = folder_dfs[0]
    for df in folder_dfs[1:]:
        folder_merged = folder_merged.merge(df, on="ID", how="outer")
    # Merge with global dataframe
    if merged_df is None:
        merged_df = folder_merged
    else:
        merged_df = merged_df.merge(folder_merged, on="ID", how="outer")

# ---------filters-----
# load filtering criteria from json
with open(args.filters, 'r') as f:
    filters = json.load(f)

# create output folders for good and bad models
good_dir = os.path.join(out_dir, "accepted")
bad_dir = os.path.join(out_dir, "failed")
os.makedirs(good_dir, exist_ok=True)
os.makedirs(bad_dir, exist_ok=True) 
accepted_metrics_file = os.path.join(good_dir, "accepted_metrics1.csv")
failed_metrics_file = os.path.join(bad_dir, "failed_metrics1.csv")


# apply filters and move files
for index, row in merged_df.iterrows():
    model_id = row['ID']
    is_good = True
    for metric, thres_cond in filters.items():
        threshold = thres_cond['threshold']
        condition = thres_cond['condition']
        if condition=="less" and metric in row and row[metric] >= threshold:
            is_good = False
            print(f"Filtering out {model_id} due to {metric}={row[metric]} >= {threshold}")
            break
        elif condition=="greater" and metric in row and row[metric] <= threshold:
            is_good = False
            print(f"Filtering out {model_id} due to {metric}={row[metric]} <= {threshold}")
            break
        else:
            print(f"Warning: Unknown condition '{condition}' for metric '{metric}' or metric not found in row.")

    # move files and data based on filtering results
    src_folder = os.path.join(out_dir, model_id)
    if is_good:
        dest_folder = os.path.join(good_dir, model_id)
        good_metrics = merged_df.iloc[[index]]
        good_metrics.to_csv(accepted_metrics_file, mode='a', header=not os.path.exists(accepted_metrics_file), index=False)
    else:
        dest_folder = os.path.join(bad_dir, model_id)
        bad_metrics = merged_df.iloc[[index]]
        bad_metrics.to_csv(failed_metrics_file, mode='a', header=not os.path.exists(accepted_metrics_file), index=False)

    if os.path.exists(src_folder):
        os.rename(src_folder, dest_folder)

# --- SAVE ---
if merged_df is not None:
    output_file = os.path.join(out_dir, "merged_af3_ter_scoring.csv")
    merged_df.to_csv(output_file, index=False)
    print(f"Metrics saved to {output_file}")
else:
    print("No CSV files found.")



