"""Training AUROC, validation-selected thresholds, full-dataset metrics.

Produces:
  results/human_train_val_metrics.csv -> Table 4 (train - test AUROC), Table S4 (validation thresholds)
  results/human_full05.csv            -> Table S1 and Table 4 ("Full-data best -> test rank")
  results/mouse_metrics.csv           -> Table 7 (train - test AUROC), Table S3 (pilot mouse sets)
"""
import sys
import warnings
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from common.data_io import load_predictions
from common.metrics import best_threshold, mets

warnings.filterwarnings("ignore")


def main():
    T = config.THRESHOLD
    hv, full = [], []
    for mod in config.HUMAN_MODS:
        for m, d in load_predictions("human", mod).items():
            y = d.label.values; s = d.score.values
            te = d.split.values == "test"; va = d.split.values == "val"; tr = d.split.values == "train"
            tv = best_threshold(y[va], s[va]); Tv = mets(y[te], s[te], tv)
            hv.append(dict(mod=mod, method=m, train_AUC=roc_auc_score(y[tr], s[tr]),
                           val_thr=tv, Tv_MCC=Tv["MCC"], Tv_Sn=Tv["Sn"], Tv_Sp=Tv["Sp"],
                           T_AP=average_precision_score(y[te], s[te])))
            F = mets(y, s, T)
            full.append(dict(mod=mod, method=m, **{k: F[k] for k in ["Acc", "Pre", "Sn", "Sp", "F1", "MCC"]},
                             AUC=roc_auc_score(y, s)))
    pd.DataFrame(hv).to_csv(config.OUT_RES / "human_train_val_metrics.csv", index=False)
    pd.DataFrame(full).to_csv(config.OUT_RES / "human_full05.csv", index=False)

    mo = []
    for mod in config.MOUSE_MODS:
        for m, d in load_predictions("mouse", mod).items():
            y = d.label.values; s = d.score.values
            te = d.split.values == "test"; tr = d.split.values == "train"
            F = mets(y, s, T); Te = mets(y[te], s[te], T)
            mo.append(dict(mod=mod, method=m, n=len(y), n_test=int(te.sum()), full_MCC=F["MCC"],
                           full_AUC=roc_auc_score(y, s), train_AUC=roc_auc_score(y[tr], s[tr]),
                           **{f"T_{k}": v for k, v in Te.items()}, T_AUC=roc_auc_score(y[te], s[te]),
                           T_AP=average_precision_score(y[te], s[te])))
    pd.DataFrame(mo).to_csv(config.OUT_RES / "mouse_metrics.csv", index=False)
    print("written: human_train_val_metrics.csv, human_full05.csv, mouse_metrics.csv")


if __name__ == "__main__":
    main()
