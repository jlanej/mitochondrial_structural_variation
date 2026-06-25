"""Pure-Python statistics kernel for the LOD analysis (no numpy/scipy).

Ports MitoHPC's lod_report.py methodology so it runs anywhere:
  * Wilson 95% score intervals for per-cell detection proportions,
  * Firth-penalized logistic dose-response fit (separation-robust — the
    recommended fix over MitoHPC's Nelder-Mead on the raw likelihood, which is
    unstable on near-separable low-replicate grids),
  * LOD50 / LOD95 by inverting the fitted curve, with a cluster (replicate-unit)
    bootstrap CI on LOD95,
  * an empirical transition read-out (the headline, robust when the fit is
    near-separable),
  * ROC / PR (AUROC trapezoid, AUPRC step) + MCC for ranking caller confidence.

Heteroplasmy is always a fraction (0..1).
"""
from __future__ import annotations

import math
import zlib

NAN = float("nan")


def expit(x):
    if x < -700:
        return 0.0
    if x > 700:
        return 1.0
    return 1.0 / (1.0 + math.exp(-x))


def logit(p):
    p = min(max(p, 1e-12), 1 - 1e-12)
    return math.log(p / (1.0 - p))


def wilson(k, n, z=1.96):
    """Point proportion + Wilson score 95% CI (clamped to [0,1])."""
    if n <= 0:
        return (NAN, NAN, NAN)
    p = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))) / denom
    return (p, max(0.0, center - half), min(1.0, center + half))


# --------------------------------------------------------------------------- #
# Firth-penalized logistic regression (2 params: intercept + slope on x=vaf)
# --------------------------------------------------------------------------- #
def _inv2(a, b, c, d):
    det = a * d - b * c
    if abs(det) < 1e-300:
        return None
    return (d / det, -b / det, -c / det, a / det)


def fit_logistic_firth(xs, ys, maxit=200, tol=1e-8):
    """Firth-penalized binary logistic fit. Returns (b0, b1) or None.

    Robust to perfect/quasi separation (returns finite, shrunk estimates) which
    a plain MLE cannot. xs = heteroplasmy fractions, ys = 0/1 detection.
    """
    n = len(xs)
    if n < 3 or len(set(xs)) < 2:
        return None
    b0, b1 = 0.0, 0.0
    for _ in range(maxit):
        # X'WX (2x2), score with Firth adjustment
        s11 = s12 = s22 = 0.0
        for x in xs:
            p = expit(b0 + b1 * x)
            w = max(p * (1 - p), 1e-12)
            s11 += w
            s12 += w * x
            s22 += w * x * x
        inv = _inv2(s11, s12, s12, s22)
        if inv is None:
            return None
        i11, i12, i21, i22 = inv
        u0 = u1 = 0.0
        for x, y in zip(xs, ys):
            p = expit(b0 + b1 * x)
            w = max(p * (1 - p), 1e-12)
            # hat value h = w * v' (X'WX)^-1 v, v=[1,x]
            v_inv_0 = i11 + i12 * x
            v_inv_1 = i21 + i22 * x
            h = w * (v_inv_0 + v_inv_1 * x)
            adj = (y - p) + h * (0.5 - p)
            u0 += adj
            u1 += adj * x
        d0 = i11 * u0 + i12 * u1
        d1 = i21 * u0 + i22 * u1
        b0 += d0
        b1 += d1
        if abs(d0) + abs(d1) < tol:
            break
    if not (math.isfinite(b0) and math.isfinite(b1)):
        return None
    return (b0, b1)


def lod_at(beta, p):
    """Heteroplasmy at which fitted P(detect) == p; NaN if slope ~0."""
    if beta is None:
        return NAN
    b0, b1 = beta
    if abs(b1) < 1e-9:
        return NAN
    return (logit(p) - b0) / b1


def cluster_bootstrap_lod(units, p=0.95, B=400, seed=7):
    """Percentile 95% CI for LOD at p by resampling replicate UNITS.

    units = list of (vaf, detected 0/1) over all levels for one (caller,variant,
    depth). Returns (lo, hi) in [0,1] or (NaN, NaN) if unstable.
    """
    n = len(units)
    if n < 8:
        return (NAN, NAN)
    rng = _Rng(seed)
    vals, fails = [], 0
    for _ in range(B):
        idx = [rng.randint(0, n - 1) for _ in range(n)]
        xs = [units[i][0] for i in idx]
        ys = [units[i][1] for i in idx]
        beta = fit_logistic_firth(xs, ys)
        if beta is None:
            fails += 1
            continue
        v = lod_at(beta, p)
        if math.isfinite(v) and 0.0 <= v <= 1.0:
            vals.append(v)
    if fails > B // 2 or len(vals) < 10:
        return (NAN, NAN)
    vals.sort()
    return (_pct(vals, 2.5), _pct(vals, 97.5))


def empirical_lod(levels):
    """levels = list of (vaf, rate) at vaf>0, sorted ascending by vaf.

    Returns dict: transition_hi (max vaf with rate<0.5), reliable_lo (min vaf
    with rate>=0.9), near_separable (<=1 partially-detected level).
    """
    pos = [(v, r) for (v, r) in sorted(levels) if v > 0]
    trans = max([v for (v, r) in pos if r < 0.5], default=NAN)
    reliable = min([v for (v, r) in pos if r >= 0.9], default=NAN)
    return {"transition_hi": trans, "reliable_lo": reliable}


def empirical_separable(counts):
    """counts = list of (k, n); near-separable if <=1 cell is 0<k<n."""
    partial = sum(1 for (k, n) in counts if 0 < k < n)
    return partial <= 1


# --------------------------------------------------------------------------- #
# ROC / PR / MCC for ranking caller confidence (when a score is available)
# --------------------------------------------------------------------------- #
def roc_pr(scores, labels):
    """scores higher = more confident; labels 1=true, 0=false. Returns dict."""
    pts = sorted(zip(scores, labels), key=lambda t: -t[0])
    P = sum(labels)
    N = len(labels) - P
    if P == 0 or N == 0:
        return None
    tp = fp = 0
    roc = [(0.0, 0.0)]
    pr = []
    for s, y in pts:
        if y == 1:
            tp += 1
        else:
            fp += 1
        roc.append((fp / N, tp / P))
        pr.append((tp / P, tp / (tp + fp)))  # (recall, precision)
    auroc = 0.0
    for i in range(1, len(roc)):
        auroc += (roc[i][0] - roc[i - 1][0]) * (roc[i][1] + roc[i - 1][1]) / 2.0
    auprc, prev_r = 0.0, 0.0
    for (r, pr_) in pr:
        auprc += (r - prev_r) * pr_
        prev_r = r
    return {"auroc": auroc, "auprc": auprc, "prevalence": P / (P + N),
            "roc": roc, "pr": pr}


def mcc(tp, fp, tn, fn):
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return 0.0 if denom == 0 else (tp * tn - fp * fn) / denom


# --------------------------------------------------------------------------- #
# small helpers (deterministic RNG + percentile) — keep stdlib-only
# --------------------------------------------------------------------------- #
class _Rng:
    """Deterministic LCG-seeded RNG for reproducible bootstrap (stdlib random
    would also work; this keeps it independent of PYTHONHASHSEED nuances)."""
    def __init__(self, seed):
        import random as _r
        self._r = _r.Random(zlib.crc32(str(seed).encode()) & 0x7FFFFFFF)

    def randint(self, a, b):
        return self._r.randint(a, b)


def _pct(sorted_vals, q):
    if not sorted_vals:
        return NAN
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = (q / 100.0) * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def boxstats(vals):
    """Five-number summary + Tukey whiskers/outliers + mean for a boxplot.

    Returns None for empty input. Whiskers extend to the most extreme value
    within 1.5*IQR of the quartiles; points beyond are listed as outliers.
    """
    xs = sorted(v for v in vals if v is not None and v == v)
    n = len(xs)
    if n == 0:
        return None
    q1, med, q3 = _pct(xs, 25), _pct(xs, 50), _pct(xs, 75)
    iqr = q3 - q1
    lo_fence, hi_fence = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    inside = [x for x in xs if lo_fence <= x <= hi_fence]
    return {
        "n": n, "min": xs[0], "q1": q1, "med": med, "q3": q3, "max": xs[-1],
        "wlo": min(inside) if inside else xs[0],
        "whi": max(inside) if inside else xs[-1],
        "mean": sum(xs) / n,
        "outliers": [x for x in xs if x < lo_fence or x > hi_fence],
    }
