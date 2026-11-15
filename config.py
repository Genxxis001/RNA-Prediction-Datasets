"""Shared settings for every script in this package.

Edit DATA_HUMAN / DATA_MOUSE if your data lives elsewhere.
Expected layout (created by prepare_data.py):
  data/human/fasta/human_hg38_<MOD>_balanced_<N>.fasta
  data/human/predictions/<MOD>/<METHOD>.csv   (id,label,score,split)
  data/human/splits/<MOD>_split.csv
  data/mouse/fasta/mouse_mm10_<MOD>_balanced_<N>.fasta
  data/mouse/predictions/<MOD>/<METHOD>.csv
  data/mouse/splits/<MOD>_split.csv
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_HUMAN = ROOT / "data" / "human"
DATA_MOUSE = ROOT / "data" / "mouse"
OUT = ROOT / "outputs"
OUT_FIG = OUT / "figures"        # PNGs, same file names as in the Overleaf project
OUT_RES = OUT / "results"        # CSV / JSON numbers behind every table and figure
OUT_TAB = OUT / "latex_tables"   # LaTeX code of Tables 3, 4, 5, 7, S1, S3, S4

# Modification order used for computations (kept identical to the original
# analysis so that the bootstrap random stream gives the same CIs).
HUMAN_MODS = ["m6A", "Nm", "m1A", "m7G", "m5C", "Pseudo"]
MOUSE_MODS = ["m6A", "m5C", "Pseudo", "m1A", "Nm", "m7G"]
MOUSE_FULL_MODS = ["m6A", "m5C", "Pseudo"]      # full-size mouse sets
MOUSE_PILOT_MODS = ["m1A", "m7G", "Nm"]         # pilot sets (Table S3)
TABLE_ORDER = ["m6A", "m1A", "m5C", "m7G", "Pseudo", "Nm"]  # order in the paper

THRESHOLD = 0.5      # common decision threshold for all predictors
N_BOOT = 1000        # bootstrap resamples for 95% CIs
SEED = 42            # random seed
WINDOW = 41          # sequence window length (nt)
BIAS_POINTS = 10     # |Sn - Sp| >= 10 points -> "biased" predictor (Fig. 3, Table 5)

# Display names
MOD_LABEL = {"Pseudo": "Ψ"}
METHOD_SHORT = {"CLSM6A-41nt": "CLSM6A", "deepSRAMP-seqonly41": "deepSRAMP"}
COLORS5 = ["#1f5fbf", "#e69500", "#2a9d3a", "#a23cc8", "#c0392b"]
BASE_COLORS = {"A": "#2a9d3a", "C": "#1f5fbf", "G": "#e69500", "U": "#a23cc8"}

for d in (OUT_FIG, OUT_RES, OUT_TAB):
    d.mkdir(parents=True, exist_ok=True)
