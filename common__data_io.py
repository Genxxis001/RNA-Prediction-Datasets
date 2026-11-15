"""Reading FASTA files and prediction CSVs, with consistency checks."""
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402


def species_dir(species):
    return config.DATA_HUMAN if species == "human" else config.DATA_MOUSE


def find_fasta(species, mod):
    tag = "hg38" if species == "human" else "mm10"
    hits = glob.glob(str(species_dir(species) / "fasta" / f"*{tag}_{mod}_balanced_*.fasta"))
    if not hits:
        raise FileNotFoundError(f"No FASTA for {species} {mod}. Run prepare_data.py first.")
    return hits[0]


def read_fasta(path):
    """Return (headers, sequences as RNA, labels) from a FASTA whose headers end with label=0/1."""
    H, S, L = [], [], []
    h = None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            h = line[1:]
            continue
        if line:
            H.append(h)
            S.append(line.upper().replace("T", "U"))
            L.append(int(h.split("label=")[1].split("|")[0]))
    return H, S, np.array(L)


def method_files(species, mod):
    return sorted(glob.glob(str(species_dir(species) / "predictions" / mod / "*.csv")))


def load_predictions(species, mod, check=True):
    """Dict method -> DataFrame(id,label,score,split). Checks ids/labels/splits against FASTA."""
    D = {os.path.basename(f)[:-4]: pd.read_csv(f) for f in method_files(species, mod)}
    if check:
        hdr, _, lab = read_fasta(find_fasta(species, mod))
        sp = pd.read_csv(species_dir(species) / "splits" / f"{mod}_split.csv")
        assert list(sp.id) == hdr, f"{species} {mod}: split ids differ from FASTA"
        for m, d in D.items():
            assert list(d.id) == hdr, f"{m}: ids differ from FASTA"
            assert (d.split.values == sp.split.values).all(), f"{m}: split differs"
            assert (d.label.values == lab).all(), f"{m}: labels differ from FASTA"
    return D
