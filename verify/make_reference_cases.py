# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Build the self-contained reference test set shipped with the package.

Synthetic, seeded images (no patient data), chosen to exercise every code path:
  structured  brain-like MRI / functional blob / weighted fusion (256x256)
  noise       uniform random images
  saturated   many 0 and 255 pixels (uint8 saturation in Q_MI)
  low_range   all values in [0, 40] (min-max normalisation in Q_NCIE)
  unit_range  fused image with values in {0, 1} (the [0,1] input branch of Q_MLI/Q_AG/Q_En)
  odd_size    255x257 (odd sizes in Q_M / dwt2, SF, gradients)
  large_512   512x512 (the down-sampling branch of Q_FSIM)
For Q_FMI the core is checked on 64x64 Sobel edge maps (Octave's edge() differs from MATLAB's),
which are stored too.

    python make_reference_cases.py        # -> reference_cases.mat (input for run_reference.m)
"""
import os
import sys

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from fusion_metrics import edge_sobel  # noqa: E402


def blob(h, w, cy, cx, s, amp):
    y, x = np.mgrid[0:h, 0:w]
    return amp * np.exp(-((y - cy) ** 2 + (x - cx) ** 2) / (2 * s * s))


def structured(h, w, rng):
    y, x = np.mgrid[0:h, 0:w]
    r = np.hypot(y - h / 2, x - w / 2)
    mri = 120 * (r < 0.42 * min(h, w)) + 60 * np.sin(r / 4) * (r < 0.38 * min(h, w)) + rng.normal(0, 4, (h, w))
    fun = blob(h, w, h * 0.4, w * 0.55, h / 9, 230) + blob(h, w, h * 0.62, w * 0.35, h / 14, 160)
    mri, fun = np.clip(mri, 0, 255), np.clip(fun + rng.normal(0, 2, (h, w)), 0, 255)
    fus = np.clip(0.6 * mri + 0.5 * fun, 0, 255)
    return [np.floor(v + 0.5).astype(np.uint8) for v in (mri, fun, fus)]


def cases():
    rng = np.random.default_rng(2026)
    C = {}
    C["structured"] = structured(256, 256, rng)
    a, b = (rng.integers(0, 256, (256, 256)).astype(np.uint8) for _ in range(2))
    C["noise"] = [a, b, ((a.astype(int) + b) // 2).astype(np.uint8)]
    m, f, _ = structured(256, 256, rng)
    C["saturated"] = [np.where(m > 100, 255, m).astype(np.uint8), np.where(f > 120, 255, 0).astype(np.uint8),
                      np.where(f > 120, 255, m).astype(np.uint8)]
    m, f, u = structured(256, 256, rng)
    C["low_range"] = [(m // 7).astype(np.uint8), (f // 7).astype(np.uint8), (u // 7).astype(np.uint8)]
    C["unit_range"] = [m, f, (u > 100).astype(np.uint8)]
    C["odd_size"] = structured(255, 257, rng)
    C["large_512"] = structured(512, 512, rng)
    return C


def main():
    C = cases()
    out = {"names": np.array(list(C), dtype=object)}
    for i, (n, (a, b, f)) in enumerate(C.items(), 1):
        out[f"A{i}"], out[f"B{i}"], out[f"F{i}"] = a, b, f
        crop = (slice(96, 160), slice(96, 160))
        for k, X in (("EA", a), ("EB", b), ("EF", f)):
            out[f"{k}{i}"] = edge_sobel(X[crop]).astype(np.uint8)
    sio.savemat(os.path.join(HERE, "reference_cases.mat"), out, do_compression=True)
    print(len(C), "cases:", ", ".join(C))


if __name__ == "__main__":
    main()
