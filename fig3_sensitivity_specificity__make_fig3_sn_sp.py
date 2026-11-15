"""Figure 3: sensitivity vs specificity of the 30 human predictors on the held-out test set
(threshold 0.5). Predictors with |Sn - Sp| >= 10 points are highlighted in blue and labelled.

Needs: outputs/results/final_test05.csv (run tables_and_statistics/compute_test_metrics.py first)
Output: outputs/figures/fig3_sn_sp.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config


def main():
    R = pd.read_csv(config.OUT_RES / "final_test05.csv"); H = R[R.sp == "human"].copy()
    H["b"] = (H.Sn - H.Sp).abs() >= config.BIAS_POINTS
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    ax.plot([45, 100], [45, 100], ls="--", color="#999", lw=0.9, zorder=0)
    g = H[~H.b]; ax.scatter(g.Sp, g.Sn, s=42, color="#c9c6bd", zorder=2)
    b = H[H.b]; ax.scatter(b.Sp, b.Sn, s=55, color="#1f5fbf", zorder=3)
    off = {"m6A-IIN": (-8, 0, "right"), "BiRNA-BERT": (-8, 0, "right"), "m5C-pred": (8, 5, "left")}
    for _, r in b.iterrows():
        dx, dy, ha = off.get(r.method, (8, 0, "left"))
        ax.annotate(f"{config.METHOD_SHORT.get(r.method, r.method)} ({config.MOD_LABEL.get(r['mod'], r['mod'])})",
                    (r.Sp, r.Sn), xytext=(dx, dy), textcoords="offset points", va="center", ha=ha, fontsize=9)
    ax.set_xlim(10, 100); ax.set_ylim(45, 100)
    ax.set_xlabel("Specificity (%)"); ax.set_ylabel("Sensitivity (%)")
    ax.grid(color="#eee", lw=0.8)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(config.OUT_FIG / "fig3_sn_sp.png", dpi=200)
    print(f"written fig3_sn_sp.png ({len(b)} biased predictors highlighted)")


if __name__ == "__main__":
    main()
