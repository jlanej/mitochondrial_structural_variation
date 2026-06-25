#!/usr/bin/env python3
"""Unit tests for the LOD statistics kernel (pipeline/lod/lod_stats.py).

Runs without Docker or scipy/numpy. Validates Wilson CIs, the Firth-penalized
logistic LOD recovery (incl. the separation case that breaks a plain MLE),
empirical LOD, and ROC/PR.

Run:  python3 test/test_lod_stats.py   (or: pytest test/test_lod_stats.py)
"""
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline", "lod"))
import lod_stats as S  # noqa: E402


def test_wilson():
    p, lo, hi = S.wilson(8, 10)
    assert abs(p - 0.8) < 1e-9 and lo < 0.8 < hi and lo > 0.4
    p, lo, hi = S.wilson(0, 10)
    assert p == 0 and lo == 0 and 0.2 < hi < 0.35    # one-sided-ish upper bound
    assert S.wilson(0, 0)[0] != S.wilson(0, 0)[0]    # NaN for n=0


def test_firth_recovers_lod():
    rng = random.Random(1)
    center = 0.05
    xs, ys = [], []
    for v in [0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.20]:
        for _ in range(30):
            xs.append(v)
            ys.append(1 if rng.random() < S.expit(60 * (v - center)) else 0)
    beta = S.fit_logistic_firth(xs, ys)
    assert beta is not None
    lod50 = S.lod_at(beta, 0.5)
    assert abs(lod50 - center) < 0.02, lod50           # within 2% of true 5%
    assert S.lod_at(beta, 0.95) > lod50                 # 95% needs higher het


def test_firth_handles_perfect_separation():
    # a clean step response: plain MLE diverges; Firth must return a finite fit
    xs = [0, 0, 0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2]
    ys = [0, 0, 0, 1, 1, 1, 1, 1, 1]
    beta = S.fit_logistic_firth(xs, ys)
    assert beta is not None and all(math.isfinite(b) for b in beta)
    assert 0.0 < S.lod_at(beta, 0.5) < 0.2


def test_fit_none_on_constant_x():
    assert S.fit_logistic_firth([0.1, 0.1, 0.1], [0, 1, 1]) is None


def test_empirical_lod():
    levels = [(0.01, 0.0), (0.02, 0.1), (0.03, 0.4), (0.05, 0.8), (0.08, 1.0), (0.10, 1.0)]
    e = S.empirical_lod(levels)
    assert abs(e["transition_hi"] - 0.03) < 1e-9    # highest het still <50%
    assert abs(e["reliable_lo"] - 0.08) < 1e-9      # lowest het at >=90%
    assert S.empirical_separable([(0, 30), (30, 30), (30, 30)]) is True
    assert S.empirical_separable([(0, 30), (10, 30), (20, 30), (30, 30)]) is False


def test_bootstrap_ci_brackets_lod():
    rng = random.Random(2)
    units = []
    for v in [0, 0.02, 0.03, 0.05, 0.08, 0.10, 0.20, 0.30]:
        for _ in range(20):
            units.append((v, 1 if (v > 0 and rng.random() < S.expit(70 * (v - 0.05))) else 0))
    lo, hi = S.cluster_bootstrap_lod(units, 0.95, B=150)
    assert math.isfinite(lo) and math.isfinite(hi) and 0 <= lo <= hi <= 1


def test_roc_pr_and_mcc():
    r = S.roc_pr([0.9, 0.8, 0.7, 0.6, 0.2, 0.1], [1, 1, 0, 1, 0, 0])
    assert r is not None and 0.5 < r["auroc"] <= 1 and 0.5 < r["auprc"] <= 1
    assert r["prevalence"] == 0.5
    assert S.roc_pr([1, 2, 3], [1, 1, 1]) is None       # one class -> None
    assert abs(S.mcc(5, 0, 5, 0) - 1.0) < 1e-9
    assert S.mcc(0, 0, 0, 0) == 0.0


def _main():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn(); print("PASS", fn.__name__)
        except AssertionError as e:
            failed += 1; print("FAIL", fn.__name__, "->", e)
        except Exception as e:  # noqa: BLE001
            failed += 1; print("ERROR", fn.__name__, "->", repr(e))
    print("\n%d/%d passed" % (len(fns) - failed, len(fns)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_main())
