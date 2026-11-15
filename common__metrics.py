"""Threshold metrics, bootstrap CI, DeLong test (Sun & Xu 2014) and Holm correction."""
import numpy as np
from scipy.stats import norm
from sklearn.metrics import matthews_corrcoef


def mets(y, s, t=0.5):
    """Acc, Pre, Sn, Sp, F1 (in %), MCC and confusion counts at threshold t (score >= t -> positive)."""
    p = (s >= t).astype(int)
    TP = int(((p == 1) & (y == 1)).sum()); TN = int(((p == 0) & (y == 0)).sum())
    FP = int(((p == 1) & (y == 0)).sum()); FN = int(((p == 0) & (y == 1)).sum())
    sn = TP / (TP + FN); sp = TN / (TN + FP); pre = TP / (TP + FP) if TP + FP else 0
    f1 = 2 * pre * sn / (pre + sn) if pre + sn else 0
    return dict(Acc=(TP + TN) / len(y) * 100, Pre=pre * 100, Sn=sn * 100, Sp=sp * 100, F1=f1 * 100,
                MCC=matthews_corrcoef(y, p), TP=TP, FN=FN, TN=TN, FP=FP)


def bootstrap_mcc_ci(y, p, rng, n_boot=1000):
    """Percentile 95% CI of MCC over n_boot resamples of the test records (fixed predictions p)."""
    bm = []
    for _ in range(n_boot):
        b = rng.integers(0, len(y), len(y))
        bm.append(matthews_corrcoef(y[b], p[b]))
    return np.percentile(bm, 2.5), np.percentile(bm, 97.5)


def best_threshold(y, s):
    """Threshold that maximises MCC (used on the validation partition for Table S4)."""
    cand = np.unique(s)
    if len(cand) > 2000:
        cand = np.quantile(s, np.linspace(0, 1, 2001))
    return float(max(cand, key=lambda t: matthews_corrcoef(y, (s >= t).astype(int))))


def _midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N); T2[J] = T
    return T2


def delong(y, s1, s2):
    """Paired DeLong test for two correlated AUROCs. Returns (AUC1 - AUC2, two-sided p)."""
    pos = s1[y == 1], s2[y == 1]; neg = s1[y == 0], s2[y == 0]
    m = len(pos[0]); n = len(neg[0]); aucs = []; V10 = []; V01 = []
    for k in range(2):
        tx = _midrank(pos[k]); ty = _midrank(neg[k]); tz = _midrank(np.concatenate([pos[k], neg[k]]))
        aucs.append((tz[:m].sum() - m * (m + 1) / 2) / (m * n))
        V10.append((tz[:m] - tx) / n); V01.append(1 - (tz[m:] - ty) / m)
    S = np.cov(np.vstack(V10)) / m + np.cov(np.vstack(V01)) / n
    z = (aucs[0] - aucs[1]) / np.sqrt(S[0, 0] + S[1, 1] - 2 * S[0, 1])
    return aucs[0] - aucs[1], 2 * norm.sf(abs(z))


def holm(items):
    """items: list of (name, diff, p). Returns list of (name, diff, p, p_holm), sorted by p."""
    items = sorted(items, key=lambda x: x[2]); k = len(items); run = 0; out = []
    for i, (m, d, p) in enumerate(items):
        run = max(run, min(1, (k - i) * p))
        out.append((m, d, p, run))
    return out
