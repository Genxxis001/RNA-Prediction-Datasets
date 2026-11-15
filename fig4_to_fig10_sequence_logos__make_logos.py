"""Figures 4-6 (human) and 8-10 (mouse): two-sample logos + top-30 ENAC features,
plus the flank-GC comparison and the numbers summarised in Table 6.

Panel A - two-sample logo: for every position (-20..+20) and base, the frequency in positives
  minus the frequency in negatives (percentage points). Only differences significant in a
  two-proportion z-test after Bonferroni correction over 41 x 4 = 164 tests (p < 0.05/164)
  are drawn; letters above 0 are enriched in positives, below 0 depleted.
Panel B - ENAC (enhanced nucleic-acid composition, window 2 -> 40 windows x 4 bases = 160
  features), ranked by information gain (mutual information with the label, bits); top 30 shown.
GC - GC content of the flanks (window minus the central 11 nt), positives vs negatives,
  Mann-Whitney U test.

Outputs: outputs/figures/<name>.png (names below), outputs/results/sequence_summary_<species>.json,
         outputs/results/flank_gc.csv
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from scipy.stats import mannwhitneyu, norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from common.data_io import find_fasta, read_fasta

B = "ACGU"; COL = config.BASE_COLORS
FP = FontProperties(family="DejaVu Sans", weight="bold")
# paper figure  ->  file name used in the Overleaf project
FILES = {
    "human": {"m6A": "fig4a_logo_m6A", "m1A": "fig4b_logo_m1A",      # Figure 4a / 4b
              "m5C": "figS1_logo_m5C", "m7G": "figS2_logo_m7G",      # Figure 5a / 5b
              "Pseudo": "figS3_logo_Psi", "Nm": "figS4_logo_Nm"},    # Figure 6a / 6b
    "mouse": {"m6A": "figS5_logo_mouse_m6A", "m5C": "figS6_logo_mouse_m5C",   # Figure 8a / 8b
              "Pseudo": "figS7_logo_mouse_Psi", "m1A": "figS8_logo_mouse_m1A",  # Figure 9a / 9b
              "m7G": "figS9_logo_mouse_m7G", "Nm": "figS10_logo_mouse_Nm"},  # Figure 10a / 10b
}


def onehot(S):
    X = np.zeros((len(S), config.WINDOW, 4))
    for i, s in enumerate(S):
        for j, c in enumerate(s):
            if c in B:
                X[i, j, B.index(c)] = 1
    return X


def letter(ax, c, x, y, w, h):
    """Draw base c as a glyph stretched to height h at (x, y)."""
    if h <= 0:
        return
    tp = TextPath((0, 0), c, size=1, prop=FP); bb = tp.get_extents()
    t = Affine2D().translate(-bb.x0, -bb.y0).scale(w * 0.9 / bb.width, h / bb.height).translate(x - w * 0.45, y)
    ax.add_patch(PathPatch(t.transform_path(tp), fc=COL[c], ec="none"))


def mi(f, y):
    """Information gain (bits) of a discrete feature f about the binary label y."""
    H = lambda p: -sum(q * math.log2(q) for q in p if q > 0)
    out = H(np.bincount(y, minlength=2) / len(y))
    for v in np.unique(f):
        m = f == v
        out -= m.mean() * H(np.bincount(y[m], minlength=2) / m.sum())
    return out


def analyse(species, mod):
    _, S, y = read_fasta(find_fasta(species, mod))
    assert all(len(s) == config.WINDOW for s in S)
    label = ("mouse " if species == "mouse" else "") + config.MOD_LABEL.get(mod, mod)
    X = onehot(S); P = X[y == 1].mean(0); N = X[y == 0].mean(0)
    n1 = int((y == 1).sum()); n0 = int((y == 0).sum())
    # two-proportion z-test per position/base, Bonferroni over 164 tests
    pp = X.sum(0) / len(y); se = np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n0)); se[se == 0] = np.inf
    z = (P - N) / se; p = 2 * norm.sf(abs(z)); sig = p < 0.05 / (config.WINDOW * 4)
    D = np.where(sig, P - N, 0) * 100
    # ENAC, window 2
    E = (X[:, :-1, :] + X[:, 1:, :]) / 2
    ig = sorted(((mi(E[:, j, b], y), j - 20, b) for j in range(E.shape[1]) for b in range(4)), reverse=True)
    top = ig[:30]
    pos = list(range(-20, 21))

    fig = plt.figure(figsize=(12, 10.5)); gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.35], hspace=0.25)
    ax = fig.add_subplot(gs[0])
    for j in range(41):
        yv = 0
        for v, c in sorted([(D[j, b], B[b]) for b in range(4) if D[j, b] > 0]):
            letter(ax, c, pos[j], yv, 0.9, v); yv += v
        yv = 0
        for v, c in sorted([(-D[j, b], B[b]) for b in range(4) if D[j, b] < 0]):
            letter(ax, c, pos[j], yv - v, 0.9, v); yv -= v
    ax.axhline(0, color="#444", lw=0.8); ax.set_xlim(-20.6, 20.6)
    lim = max(abs(D.clip(min=0).sum(1)).max(), abs(D.clip(max=0).sum(1)).max()) * 1.08
    ax.set_ylim(-lim, lim); ax.set_xticks(pos); ax.set_xticklabels([str(x) for x in pos], fontsize=7)
    ax.set_ylabel("Enriched (+) / depleted (−)\nin positives (% points)")
    ax.set_xlabel("Position relative to modified site")
    ax.set_title(f"A  Two-sample logo, {label} (n = {n1:,} positives vs {n0:,} negatives)", loc="left", fontsize=11)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    bx = fig.add_subplot(gs[1])
    for r, (v, j, b) in enumerate(top):
        bx.plot([j - 0.45, j + 1.45], [r, r], color=COL[B[b]], lw=7, solid_capstyle="butt")
        bx.text(j + 1.8, r, f"Top{r+1}: {B[b]}@[{j}, {j+1}]  IG={v:.3f}", va="center", fontsize=7.5, clip_on=False,
                zorder=3, bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.9))
    bx.set_ylim(29.8, -0.8); bx.set_xlim(-20.6, 20.6); bx.set_yticks([])
    bx.axvline(0, color="#bbb", lw=0.8, ls="--", zorder=0)
    bx.set_xticks(pos); bx.set_xticklabels([str(x) for x in pos], fontsize=7)
    bx.set_xlabel("Position relative to modified site (ENAC window, size 2)")
    bx.set_title(f"B  Top 30 ENAC features ranked by information gain, {label}", loc="left", fontsize=11)
    for s in ["top", "right", "left"]:
        bx.spines[s].set_visible(False)
    fig.savefig(config.OUT_FIG / f"{FILES[species][mod]}.png", dpi=170, bbox_inches="tight"); plt.close(fig)

    # flank GC (positions outside the central 11 nt)
    gc = lambda s: sum(c in "GC" for c in s[:15] + s[26:]) / 30 * 100
    gp = [gc(s) for s, l in zip(S, y) if l == 1]; gn = [gc(s) for s, l in zip(S, y) if l == 0]
    gcrow = dict(species=species, mod=mod, gc_pos=np.mean(gp), gc_neg=np.mean(gn),
                 gc_diff=np.mean(gp) - np.mean(gn), mannwhitney_p=mannwhitneyu(gp, gn).pvalue)

    enr = sorted((round(float(D[j, b]), 1), pos[j], B[b]) for j in range(41) for b in range(4) if j != 20)
    summ = dict(n=int(len(y)), n_sig_cells=int(sig.sum()), max_abs_diff=round(float(abs(D).max()), 1),
                top30_within_pm5=sum(1 for v, j, b in top if -5 <= j and j + 1 <= 5),
                top30_base_counts={c: sum(1 for v, j, b in top if B[b] == c) for c in B},
                top5=[f"{B[b]}@[{j},{j+1}] {v:.3f}" for v, j, b in top[:5]],
                most_enriched=enr[-6:][::-1], most_depleted=enr[:6])
    return summ, gcrow


def main():
    gcrows = []
    for species in ["human", "mouse"]:
        summary = {}
        for mod in config.TABLE_ORDER:
            summary[mod], g = analyse(species, mod); gcrows.append(g)
            print(f"{species:5s} {mod:6s} -> {FILES[species][mod]}.png   flank GC diff {g['gc_diff']:+.1f}")
        json.dump(summary, open(config.OUT_RES / f"sequence_summary_{species}.json", "w"), indent=1, ensure_ascii=False)
    pd.DataFrame(gcrows).to_csv(config.OUT_RES / "flank_gc.csv", index=False)


if __name__ == "__main__":
    main()
