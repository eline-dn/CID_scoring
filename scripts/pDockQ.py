#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import tempfile
from pathlib import Path
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------
def merge_chains(model, chains_to_merge):
    for chain in chains_to_merge[1:]:
        for i, res in enumerate(model[chain]):
            res.id = (chain, res.id[1], res.id[2])
            model[chains_to_merge[0]].add(res)
        model.detach_child(chain)
    model[chains_to_merge[0]].id = "".join(chains_to_merge)
    return model
"""
usage: 
model = load_PDB("2j5l.cif.gz")
native = load_PDB("2j4w.cif.gz")

model = merge_chains(model, ["B", "C"])
native = merge_chains(native, ["L", "H"]) # automatic mapping will not work here, user must merge in the right order

# native:model chain map dictionary for two interfaces
chain_map = {"D":"A", "LH":"BC"}
# returns a dictionary containing the results and the total DockQ score
run_on_all_native_interfaces(model, native, chain_map=chain_map)
"""


def run_dockq(model, native, chain_map):
    res, total = run_on_all_native_interfaces(
    model, native, chain_map=chain_map )
    res, total = run_on_all_native_interfaces(
        model, native, chain_map=chain_map)
    return total


# ------------------------------------------------------------
# Core logic
# ------------------------------------------------------------

def score_binary(native_pdb, model_pdb):
    native = load_PDB(native_pdb)
    model  = load_PDB(model_pdb)
    return run_dockq(model, native, {"A": "A", "B": "B"})


def score_ternary(native_pdb, model_pdb):
    native = load_PDB(native_pdb, small_molecule=True)
    model  = load_PDB(model_pdb, small_molecule=True)

    scores = {}

    # A–L
    #scores["ter_dockq_A_L"] = run_dockq(
    #model, native, {"A": "A", "L": "L"})

    # B–L
    #scores["ter_dockq_B_L"] = run_dockq(
    #model, native, {"B": "B", "L": "L"})

    # A-B
    scores["ter_dockq_A_B"] = run_dockq(
        model, native, {"A": "A","B": "B"}
    )
    # (A+L)–B
    native_AL = merge_chains(native, ["A", "L"])
    model_AL  = merge_chains(model,  ["A", "L"])

    scores["ter_dockq_AL_B"] = run_dockq(
        model_AL, native_AL, {"AL": "AL", "B": "B"}
    )

    return scores


# ------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Compute DockQ for binary or ternary AF3 repredictions"
    )
    ap.add_argument(
        "--native-pdbs",
        required=True,
        nargs="+",
        help="List of reference/native PDB paths"
    )
    ap.add_argument(
        "--model-pdbs",
        required=True,
        nargs="+",
        help="List of AF3 repredicted PDB paths (same filenames)"
    )
    ap.add_argument(
        "--out-csv",
        required=True,
        help="Output CSV path"
    )
    ap.add_argument(
        "--ternary",
        action="store_true",
        help="Score ternary complexes (A,B,L). Default: binary A–B."
    )
    return ap.parse_args()


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    args = parse_args()

    native_map = {Path(p).stem: p for p in args.native_pdbs}
    #model_map  = {Path(p).stem: p for p in args.model_pdbs}
    model_map = {}
    for p in args.model_pdbs:
        stem = Path(p).stem
        """if not stem.endswith("_model"):
            raise ValueError(f"Model filename does not end with '_model': {p}")"""
        binder_id = stem.replace("_model","")  # remove "_model"
        model_map[binder_id] = p

    common_ids = sorted(native_map.keys() & model_map.keys())
    if not common_ids:
        raise SystemExit("No matching binder IDs between native and model PDBs")

    rows = []

    for bid in common_ids:
        row = {"id": bid}

        native_pdb = native_map[bid]
        model_pdb  = model_map[bid]

        if args.ternary:
            scores = score_ternary(native_pdb, model_pdb)
            row.update(scores)
        else:
            row["bin_dockq_A_B"] = score_binary(native_pdb, model_pdb)

        rows.append(row)

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    df.to_csv(args.out_csv, index=False)


if __name__ == "__main__":
    main()
