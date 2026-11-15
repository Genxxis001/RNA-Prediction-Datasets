"""Held-out test metrics at threshold 0.5, bootstrap MCC CIs, and DeLong + Holm tests.

Produces:
  results/final_test05.csv  -> Table 3 (human), Table 7 (mouse), Table 5 (bias), Figure 3
  results/delong.csv        -> Table 4 column "Significantly below best AUROC"
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from common.data_io import load_predictions
from common.metrics import bootstrap_mcc_ci, delong, holm, mets

warnings.filterwarnings("ignore")


def main():
    rng = np.random.default_rng(config.SEED)   # one stream, same order as the original analysis
    rows, dl = [], []
    for sp, mods in [("human", config.HUMAN_MODS), ("mouse", config.MOUSE_FULL_MODS)]:
        for mod in mods:
            D = load_predictions(sp, mod)
            for m, d in D.items():
                te = d.split == "test"; y = d.label.values[te]; s = d.score.values[te]
                M = mets(y, s, config.THRESHOLD)
                lo, hi = bootstrap_mcc_ci(y, (s >= config.THRESHOLD).astype(int), rng, config.N_BOOT)
                rows.append(dict(sp=sp, mod=mod, method=m, **M, MCC_lo=lo, MCC_hi=hi, AUC=roc_auc_score(y, s),
                                 full=mets(d.label.values, d.score.values, config.THRESHOLD)["MCC"]))
            te_of = lambda m: D[m].split == "test"
            aucs = {m: roc_auc_score(d.label[te_of(m)], d.score[te_of(m)]) for m, d in D.items()}
            best = max(aucs, key=aucs.get)
            y = D[best].label[te_of(best)].values; sb = D[best].score[te_of(best)].values
            tests = [(m, *delong(y, sb, D[m].score[te_of(m)].values)) for m in D if m != best]
            for m, diff, p, pa in holm(tests):
                dl.append(dict(sp=sp, mod=mod, best=best, other=m, dAUC=diff, p=p, p_holm=pa))
    R = pd.DataFrame(rows); DL = pd.DataFrame(dl)
    R.to_csv(config.OUT_RES / "final_test05.csv", index=False)
    DL.to_csv(config.OUT_RES / "delong.csv", index=False)
    print(f"final_test05.csv: {len(R)} rows   delong.csv: {len(DL)} rows")


if __name__ == "__main__":
    main()
