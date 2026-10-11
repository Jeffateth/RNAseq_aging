"""Rebuild metadata/samples.csv (one row per sequencing run) from the two raw NCBI tables.

Inputs (as fetched from NCBI, unmodified):
  metadata/runinfo_raw.csv    SRA run table (read counts, library layout, instrument, GEO sample id)
  metadata/biosample_raw.csv  BioSample records (sex, strain, mouse id, treatment) for the same runs
  metadata/srr_list.txt       the 27 assigned runs
Run:  python scripts/build_metadata.py          (writes metadata/samples.csv)
"""
from pathlib import Path

import pandas as pd

META = Path(__file__).resolve().parents[1] / "metadata"
GROUP_ORDER = ["aging 3mo", "aging 18mo", "aging 21mo", "AL", "aDR", "PBS", "YMP"]   # row order: by group, then mouse id
AGING_PROJECT = "PRJNA875066"   # aging cohort; the other BioProject (PRJNA946653) is the rejuvenation cohort


def build() -> pd.DataFrame:
    runs = pd.read_csv(META / "runinfo_raw.csv")
    bio = pd.read_csv(META / "biosample_raw.csv")
    d = runs.merge(bio, on="BioSample", suffixes=("", "_bio"))
    aging = d["BioProject"] == AGING_PROJECT
    age = d["title"].str.extract(r"(\d+)mo")[0].astype(int)            # age is only in the title, e.g. "Cerebellum Male 18mo 31ID"
    mouse = d["mouse id number"].fillna(d["mouse id_number"]).astype(int)
    condition = d["treatment"].where(~aging, "untreated aging cohort")
    out = pd.DataFrame({
        "SRR": d["Run"], "GSM": d["LibraryName"], "BioProject": d["BioProject"], "BioSample": d["BioSample"],
        "age_months": age, "sex": d["sex"], "tissue": d["tissue"],
        "experiment_group": d["experiment group"].where(~aging, "Aging"), "condition": condition,
        "mouse_id": mouse, "strain": d["strain"], "platform": d["Model"], "layout": d["LibraryLayout"],
        "n_reads_spots": d["spots"], "read_length_per_spot": d["avgLength"], "sra_size_MB": d["size_MB"], "title": d["title"],
        "group": condition.where(~aging, "aging " + age.astype(str) + "mo"),   # the 7 analysis groups
    })
    assigned = (META / "srr_list.txt").read_text().split()
    out = out[out["SRR"].isin(assigned)]
    key = out["group"].map({g: i for i, g in enumerate(GROUP_ORDER)})
    return out.assign(_k=key).sort_values(["_k", "mouse_id"]).drop(columns="_k").reset_index(drop=True)


if __name__ == "__main__":
    build().to_csv(META / "samples.csv", index=False)
    print("wrote", META / "samples.csv")
