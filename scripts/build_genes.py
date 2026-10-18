"""Rebuild metadata/genes_M39.csv (gene id -> gene symbol) from the GENCODE M39 transcript sequences.

Each transcript name in the FASTA is  ENST|ENSG|OTT gene|OTT transcript|transcript name|gene symbol|length|type|  (fields 2 and 6 are used).
Run:  python scripts/build_genes.py                      (downloads the FASTA from GENCODE, ~135 MB)
      python scripts/build_genes.py path/to/gencode.vM39.transcripts.fa.gz
"""
import gzip
import sys
import urllib.request
from pathlib import Path

import pandas as pd

URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.transcripts.fa.gz"
OUT = Path(__file__).resolve().parents[1] / "metadata/genes_M39.csv"


def build(source: str = URL) -> pd.DataFrame:
    handle = gzip.open(source if not source.startswith("http") else urllib.request.urlopen(source), "rt")
    genes = {}
    with handle as f:
        for line in f:
            if line.startswith(">"):
                p = line[1:].split("|")
                genes[p[1]] = p[5]
    return pd.DataFrame({"gene_id": genes.keys(), "symbol": genes.values()}).sort_values("gene_id")


if __name__ == "__main__":
    build(*sys.argv[1:2]).to_csv(OUT, index=False)
    print(f"wrote {OUT}")
