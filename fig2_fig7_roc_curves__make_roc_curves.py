"""Figure 2 (human, 2x3 panels) and Figure 7 (mouse, 1x3 panels): held-out test ROC curves.

For every modification the five predictors are drawn in decreasing test-AUROC order;
the legend shows the test AUROC. Only records with split == "test" are used.

Outputs: outputs/figures/fig2_roc_human.png, outputs/figures/fig5_roc_mouse.png
(file names as in the Overleaf project; fig5_roc_mouse.png is Figure 7 in the paper)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from common.data_io import load_predictions


def panel(ax, species, mod, title_prefix=""):
    res = []
    for m, d in load_predictions(species, mod).items():
        d = d[d.split == "test"]
        res.append((roc_auc_score(d.label, d.score), m, d))
    res.sort(reverse=True)
    for c, (auc, m, d) in zip(config.COLORS5, res):
        fpr, tpr, _ = roc_curve(d.label, d.score)
        ax.plot(fpr, tpr, color=c, lw=1.6, label=f"{config.METHOD_SHORT.get(m, m)} ({auc:.3f})")
    ax.plot([0, 1], [0, 1], ls="--", color="#999", lw=0.8)
    ax.set_title(f"{title_prefix}{config.MOD_LABEL.get(mod, mod)}  (test n = {len(res[0][2]):,})", loc="left", fontsize=11)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01)
    ax.set_xlabel("1 − Specificity"); ax.set_ylabel("Sensitivity")
    ax.legend(loc="lower right", fontsize=8, title="Method (AUROC)", title_fontsize=8, frameon=False)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)


def main():
    fig, axs = plt.subplots(2, 3, figsize=(13, 8.6))
    for ax, mod in zip(axs.flat, config.HUMAN_MODS):
        panel(ax, "human", mod)
    fig.tight_layout(); fig.savefig(config.OUT_FIG / "fig2_roc_human.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(13, 4.6))
    for ax, mod in zip(axs, config.MOUSE_FULL_MODS):
        panel(ax, "mouse", mod, "Mouse ")
    fig.tight_layout(); fig.savefig(config.OUT_FIG / "fig5_roc_mouse.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("written fig2_roc_human.png, fig5_roc_mouse.png")


if __name__ == "__main__":
    main()
