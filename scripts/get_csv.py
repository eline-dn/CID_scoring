# in each folder: 
# get csvs as pd df
# clean them up: 
# - ensure col "ID" from binder_id, id,... eg. 
# - ensure t2_2_7_3_t0.1_ltp0.2_dcut15.0_model0 or t2_2_7_4_t0.3_ltp0.3_dcut8.0_model becomes =>  t2_2_7_3_t0.1_ltp0.2_dcut15.0
# - put prefix on col names to ensure the source folder: 
# merge all on "ID" col (the csvs from the three folders)
# save metrics.csv in output 


from pathlib import Path
import pandas as pd
import re

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--base_folder', type=str, required=True, help=' Target+ ligand complex output folder')
args = parser.parse_args()


# --- CONFIG ---
# Path to the top-level output folder containing subfolders
base_folder = Path(args.base_folder)
# Name of the final merged CSV
output_file = base_folder / "metrics.csv"

# --- HELPER FUNCTIONS ---
def clean_id(id_str):
    """
    Clean structure ID by removing the '_model' suffix and trailing digits.
    Example:
        t2_2_7_3_t0.1_ltp0.2_dcut15.0_model0 -> t2_2_7_3_t0.1_ltp0.2_dcut15.0
    """
    return re.sub(r'_model\d*$', '', id_str)

def clean_dataframe(df, folder_prefix):
    """
    Ensure 'ID' column exists and clean it.
    Prefix other columns with folder name to avoid collisions.
    """
    # Identify the ID column: common possibilities
    id_col_candidates = [c for c in df.columns if c.lower() in ["id", "binder_id", "ID"]]
    if not id_col_candidates:
        raise ValueError("No ID column found in dataframe")
    id_col = id_col_candidates[0]
    # Standardize the ID column
    df = df.rename(columns={id_col: "ID"})
    df["ID"] = df["ID"].astype(str).apply(clean_id)
    # Prefix other columns
    other_cols = [c for c in df.columns if c != "ID"]
    df = df.rename(columns={c: f"{folder_prefix}_{c}" for c in other_cols})
    return df
def clean(df: pd.DataFrame, folder_prefix) -> pd.DataFrame:
    """
    Collapse rows whose ID ends with '_modelX' (X = single digit)
    by averaging numeric columns.
    clean id
    add prefix
    """
    df = df.copy()
    # Identify the ID column: common possibilities
    id_col_candidates = [c for c in df.columns if c.lower() in ["id", "binder_id", "ID"]]
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


# --- MAIN ---
merged_df = None

for folder in base_folder.iterdir():
    if folder.is_dir():
        folder_prefix = folder.name
        # Read all CSVs in the folder
        folder_dfs = []
        for csv_file in folder.glob("*.csv"):
            df = pd.read_csv(csv_file)
            #df = clean_dataframe(df, folder_prefix)
            df = clean(df, folder_prefix)
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

# --- SAVE ---
if merged_df is not None:
    merged_df.to_csv(output_file, index=False)
    print(f"Metrics saved to {output_file}")
else:
    print("No CSV files found.")
