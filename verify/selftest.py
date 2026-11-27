# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Self-test: Python port vs the reference values produced by the original .m files.

    python selftest.py            # exit code 0 = all metrics agree

Inputs: reference_cases.mat (7 synthetic cases, see make_reference_cases.py).
Reference: reference_values.json (written by run_reference.m in GNU Octave + this script's
--export step). A metric passes when |python - reference| <= 1e-9 * max(1, |reference|).
"""
import json
import os
import sys

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import fusion_metrics as M  # noqa: E402

TOL = 1e-9
# Ill-conditioned (case, metric) pairs: the synthetic 'saturated' case has perfectly flat 0/255
# regions where the local variances inside Q_VIFF are pure round-off (~1e-11) and sit at the
# .m file's 1e-10 thresholds; the result then depends on the summation order of filter2, which
# differs between MATLAB, Octave and SciPy. On real images Q_VIFF agrees to 1e-14.
RELAXED = {("saturated", "Q_VIFF"): 1e-4,
           # Q_CB divides two Gaussian-filtered versions of the CSF band-passed image; where the
           # denominator crosses zero (large signal-free areas) this ratio amplifies the FFT
           # round-off (~1e-16) by up to ~1e13. Fed with Octave's own FFT output the port
           # reproduces Octave to 1e-15, so the remaining differences come from the FFT library
           # (NumPy/pocketfft vs FFTW; MATLAB differs again) and are a property of the metric.
           ("saturated", "Q_CB"): 5e-2, ("unit_range", "Q_CB"): 1e-3, ("low_range", "Q_CB"): 1e-7}


def _names(arr):
    """Case / metric names from a MATLAB cell array loaded by scipy (nested object arrays)."""
    out = []
    for x in np.ravel(arr):
        while isinstance(x, np.ndarray):
            x = x.ravel()[0]
        out.append(str(x).strip())
    return out


def python_values(S, names):
    K = len(_names(S["names"]))
    out = {}
    for k in range(1, K + 1):
        A, B, F = S[f"A{k}"], S[f"B{k}"], S[f"F{k}"]
        case = _names(S["names"])[k - 1]
        row = {}
        for n in names:
            if n == "Q_FMI_core":
                row[n] = M.Q_FMI(S[f"EA{k}"].astype(float), S[f"EB{k}"].astype(float),
                                 S[f"EF{k}"].astype(float), feature="none")
            elif n.endswith("_double") or n.endswith("_unit"):     # same metric, other input class
                base, scale = (n[:-7], 1.0) if n.endswith("_double") else (n[:-5], 255.0)
                fn, src = M.METRICS[base]
                a, b, f = (X.astype(np.float64) / scale for X in (A, B, F))
                row[n] = fn(a, b, f) if src else fn(f)
            else:
                fn, src = M.METRICS[n]
                row[n] = fn(A, B, F) if src else fn(F)
        # edge maps stored in the reference set must be reproduced bit-exactly
        crop = (slice(96, 160), slice(96, 160))
        row["edge_maps_identical"] = all(np.array_equal(M.edge_sobel(X[crop]).astype(np.uint8), S[f"{e}{k}"])
                                         for e, X in (("EA", A), ("EB", B), ("EF", F)))
        out[case] = row
    return out


def main():
    S = sio.loadmat(os.path.join(HERE, "reference_cases.mat"))
    if len(sys.argv) > 1 and sys.argv[1] == "--export":        # after run_reference.m
        O = sio.loadmat(os.path.join(HERE, "reference_values.mat"))
        names = _names(O["names"])
        cases = _names(S["names"])
        ref = {c: {n: float(O["R"][i, j]) for j, n in enumerate(names)} for i, c in enumerate(cases)}
        json.dump({"source": "original .m files run in GNU Octave 10.3 (image 2.20, signal 1.4.8)",
                   "values": ref}, open(os.path.join(HERE, "reference_values.json"), "w"), indent=1)
        print("wrote reference_values.json")
        return 0
    ref = json.load(open(os.path.join(HERE, "reference_values.json")))["values"]
    names = list(next(iter(ref.values())))
    py = python_values(S, names)
    worst, ok = {}, True
    for n in names:
        errs = []
        for c in ref:
            r, p = ref[c][n], py[c][n]
            if np.isnan(r) and np.isnan(p):
                errs.append(0.0); continue
            e = abs(p - r) / max(1.0, abs(r))
            if (c, n) in RELAXED:
                relaxed_ok = e <= RELAXED[(c, n)]
                print(f"  note: {c}/{n} rel. diff {e:.1e} (ill-conditioned case, tolerance {RELAXED[(c, n)]:.0e}:"
                      f" {'OK' if relaxed_ok else 'FAIL'})")
                ok &= relaxed_ok
                continue
            errs.append(e)
        worst[n] = max(errs)
        ok &= worst[n] <= TOL
    edges = all(py[c]["edge_maps_identical"] for c in py)
    ok &= edges
    print(f"{len(ref)} cases: {', '.join(ref)}")
    for n in names:
        print(f"  {n:12s} max rel. diff {worst[n]:.1e}  {'OK' if worst[n] <= TOL else 'FAIL'}")
    print(f"  edge_sobel reproduces the stored edge maps: {'OK' if edges else 'FAIL'}")
    print("ALL OK" if ok else "SOME FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
