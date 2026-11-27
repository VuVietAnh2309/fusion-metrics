# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Regenerate docs/API.md. Every example value is computed by the installed package.

    python docs/make_api.py
"""
import inspect
import os

import numpy as np

import fusion_metrics as fm

HERE = os.path.dirname(os.path.abspath(__file__))

# name -> (what it measures, range / direction, original reference)
INFO = {
    "Q_MLI": ("Mean local intensity of the fused image.", "[0, 1]", "—"),
    "Q_En": ("Shannon entropy of the fused image (information content).", "bits, higher = more information",
             "Roberts et al., J. Appl. Remote Sens. 2008"),
    "SF": ("Spatial frequency over rows, columns and both diagonals (activity / sharpness).", "≥ 0, higher = sharper",
           "Eskicioglu & Fisher, IEEE Trans. Commun. 1995"),
    "Q_AG": ("Average gradient with forward differences on the [0, 1] scale.", "≥ 0, higher = sharper", "—"),
    "AG_VIFB": ("Average gradient as defined in the VIFB benchmark (central differences, 0–255 scale).",
                "≥ 0, higher = sharper", "Zhang et al., VIFB, CVPR Workshops 2020"),
    "Q_MI": ("Mutual information MI(A,F) + MI(B,F).", "bits, higher = better", "Qu et al., Electron. Lett. 2002"),
    "Q_NCIE": ("Nonlinear correlation information entropy between A, B and F.", "[0, 1], higher = better",
               "Wang, Shen & Jin, 2008"),
    "Q_TE_Nava": ("Normalised Tsallis-entropy mutual information, q = 0.43137.", "higher = better",
                  "Nava et al., Proc. SPIE 2007"),
    "Q_FMI": ("Feature mutual information on Sobel edge maps, window 3.", "[0, 1], higher = better",
              "Haghighat et al., Comput. Electr. Eng. 2011"),
    "Q_ABF": ("Edge-strength-and-orientation preservation Q^AB/F.", "[0, 1], higher = better",
              "Xydeas & Petrović, Electron. Lett. 2000"),
    "Q_M": ("Multiscale edge preservation (Haar wavelet, 2 levels).", "[0, 1], higher = better",
            "Wang & Liu, ICSP 2008"),
    "Q_SF": ("Relative spatial-frequency error of F against the sources.", "≈ 0 is best; sign shows over/under-sharpening",
             "Zheng et al., Inf. Fusion 2007"),
    "Q_FSIM": ("Feature similarity (FSIM(A,F) + FSIM(B,F)) / 2.", "[0, 1], higher = better",
               "Zhang et al., IEEE Trans. Image Process. 2011"),
    "Q_PSNR": ("Mean of PSNR(F, A) and PSNR(F, B).", "dB, higher = better", "—"),
    "Q_CC": ("Mean correlation coefficient between F and each source.", "[-1, 1], higher = better", "—"),
    "Q_SCD": ("Sum of the correlations of differences.", "≤ 2, higher = better",
              "Aslantas & Bendes, AEU Int. J. Electron. Commun. 2015"),
    "VIF_Sheikh": ("Visual information fidelity VIF(A,F) + VIF(B,F), pixel-domain multiscale.",
                   "≥ 0, higher = better", "Sheikh & Bovik, IEEE Trans. Image Process. 2006"),
    "Q_VIFF": ("Visual information fidelity for fusion (VIFF).", "[0, 1], higher = better",
               "Han et al., Inf. Fusion 2013"),
    "Q_CI": ("Contrast index: population standard deviation (SD) of F on the [0, 1] scale.",
             "≥ 0, higher = more contrast", "Rao, Int. J. Remote Sens. 1997"),
    "Q_EI": ("Edge intensity sqrt(sum(Sx² + Sy²)) / (M·N) with Sobel kernels, on the input scale.",
             "≥ 0, higher = stronger edges", "Rajalingam & Priya, IJESI 2018"),
    "Q_CE": ("Cross entropy 0.5·[KL(h_A‖h_F) + KL(h_B‖h_F)] of the 256-bin histograms.", "≥ 0, lower = better",
             "Bulanon et al., Biosyst. Eng. 2009"),
    "Q_NMI": ("Normalized mutual information 2·[MI(A,F)/(H(A)+H(F)) + MI(B,F)/(H(B)+H(F))].",
              "[0, 2], higher = better", "Hossny et al., Electron. Lett. 2008"),
    "MS_SSIM": ("SSIM (SSIM(A,F) + SSIM(B,F)) / 2, single-scale, 11×11 Gaussian window.", "[-1, 1], higher = better",
                "Wang et al., IEEE Trans. Image Process. 2004"),
    "Q_UIQI": ("Universal image quality index (UIQI(A,F) + UIQI(B,F)) / 2, 8×8 window. uint8 input follows MATLAB's "
               "uint8 arithmetic inside UIQI.m; pass float64 for the textbook index.", "[-1, 1], higher = better",
               "Wang & Bovik, IEEE Signal Process. Lett. 2002"),
    "Q_TE": ("Tsallis-entropy mutual information I_q(A,F) + I_q(B,F), q = 1.85.", "higher = better",
             "Cvejic, Canagarajah & Bull, Electron. Lett. 2006"),
    "Q_NABF": ("Fusion artifacts N^AB/F: weighted share of edges that F strengthens without preserving them "
               "(modified definition).", "[0, 1], lower = better", "Shreyamsha Kumar, Signal Image Video Process. 2013"),
    "Q_P": ("Phase-congruency metric: product of the masked correlations of phase congruency and of the maximum "
            "and minimum moments.", "≤ 1, higher = better", "Zhao, Laganière & Liu, IJICIC 2007"),
    "Q_S": ("Piella's fusion quality index Q_S (local SSIM weighted by the source variances).", "[-1, 1], higher = better",
            "Piella & Heijmans, ICIP 2003"),
    "Q_Sw": ("Piella's weighted fusion quality index Q_W (saliency-weighted Q_S).", "[-1, 1], higher = better",
             "Piella & Heijmans, ICIP 2003"),
    "Q_Se": ("Piella's edge-dependent fusion quality index Q_E = Q_W(edge maps)^2 (alpha = 1).",
             "[-1, 1], higher = better", "Piella & Heijmans, ICIP 2003"),
    "Q_C": ("Cvejic's metric: SSIM maps weighted by the local source–fused covariances.", "[0, 1], higher = better",
            "Cvejic et al., Int. J. Signal Process. 2005"),
    "Q_Y": ("Yang's SSIM-based metric (weighted SSIM in redundant regions, maximum SSIM elsewhere).",
            "[0, 1], higher = better", "Yang et al., Inf. Fusion 2008"),
    "Q_H": ("SSIM of 7×7 local-entropy feature maps, weighted by the sources' relative entropy; returns mean(1 − E).",
            "higher = better", "Hossny, Nahavandi & Creighton, ICIRA 2008"),
    "Q_CV": ("Chen–Varshney metric: CSF-filtered difference energy in 16×16 regions weighted by edge energy.",
             "≥ 0, lower = better", "Chen & Varshney, Inf. Fusion 2007"),
    "Q_CB": ("Chen–Blum metric: preservation of CSF-filtered local contrast.", "[0, 1], higher = better",
             "Chen & Blum, Image Vis. Comput. 2009"),
}


def example_images():
    rng = np.random.default_rng(0)
    y, x = np.mgrid[0:128, 0:128]
    a = (128 + 100 * np.sin(x / 9.0) * np.cos(y / 13.0)).clip(0, 255).astype(np.uint8)   # "anatomical"
    b = ((np.hypot(x - 64, y - 64) < 40) * 180 + rng.integers(0, 40, a.shape)).clip(0, 255).astype(np.uint8)
    f = ((a.astype(np.float64) + b) / 2).round().astype(np.uint8)                          # naive average fusion
    return a, b, f


def main():
    A, B, F = example_images()
    vals = fm.evaluate(A, B, F)
    out = ["# API reference", "",
           f"*Generated by `docs/make_api.py` for fusion-metrics {fm.__version__}. Example values are computed by the "
           "package itself.*", "",
           "## Conventions", "",
           "- `A`, `B`: source images; `F`: fused image. 2-D grayscale arrays, `uint8` or `float64` in [0, 255].",
           "- Every metric returns a Python `float`.",
           "- Read files with `fm.read_gray(path)` (MATLAB `imread` + `rgb2gray`) to reproduce MATLAB values; "
           "`fm.read_gray_cv(path)` uses OpenCV's gray conversion instead.", "",
           "## Example inputs", "",
           "All examples below use the same deterministic 128×128 images:", "",
           "```python",
           "import numpy as np, fusion_metrics as fm",
           "rng = np.random.default_rng(0)",
           "y, x = np.mgrid[0:128, 0:128]",
           "A = (128 + 100*np.sin(x/9.0)*np.cos(y/13.0)).clip(0, 255).astype(np.uint8)",
           "B = ((np.hypot(x-64, y-64) < 40)*180 + rng.integers(0, 40, A.shape)).clip(0, 255).astype(np.uint8)",
           "F = ((A.astype(np.float64) + B) / 2).round().astype(np.uint8)   # naive average fusion",
           "```", "",
           "## High-level API", "",
           "### `evaluate(A, B, F, names=None) -> dict[str, float]`", "",
           "Compute several metrics at once. `names` is a list of metric names, a preset name "
           f"({', '.join(f'`{k}`' for k in fm.PRESETS)}), or `None` for all metrics.", "",
           "```python", ">>> fm.evaluate(A, B, F, \"feature\")",
           "{" + ", ".join(f"'{k}': {vals[k]:.6f}" for k in fm.PRESETS["feature"]) + "}", "```", "",
           "### `evaluate_files(path_a, path_b, path_f, names=None, reader=None) -> dict[str, float]`", "",
           "Same as `evaluate`, reading the three image files with `reader` (default `read_gray`).", "",
           "### `read_gray(path) -> np.ndarray[uint8]` · `read_gray_cv(path)` · `rgb2gray(rgb)`", "",
           "Image readers and MATLAB's exact `rgb2gray` (coefficients 0.298936, 0.587043, 0.114021, rounded).", "",
           "### `METRICS`, `PRESETS`", "",
           "`METRICS` maps each name to `(function, needs_sources)`; `PRESETS` maps a group name to a list of metric "
           "names.", "",
           "## Metrics", ""]
    for name, (fn, src) in fm.METRICS.items():
        what, rng_, ref = INFO[name]
        sig = "(A, B, F)" if src else "(F)"
        call = f"fm.{name}(A, B, F)" if src else f"fm.{name}(F)"
        groups = [k for k, v in fm.PRESETS.items() if k != "all" and name in v]
        out += [f"### `{name}{sig} -> float`", "",
                what, "",
                f"- **Range / direction:** {rng_}",
                f"- **Reference:** {ref}",
                f"- **Presets:** {', '.join(f'`{g}`' for g in groups) or '—'}",
                f"- **Source:** `fusion_metrics/{inspect.getmodule(fn).__name__.split('.')[-1]}.py`", "",
                "```python", f">>> {call}", f"{vals[name]:.6f}", "```", ""]
    out += ["## Command line", "", "```text", "fusion-metrics A B F [--preset NAME | --metrics NAME ...] "
            "[--reader matlab|cv] [--json]", "fusion-metrics --dirs DIR_A DIR_B DIR_F [--csv FILE] ...",
            "fusion-metrics --list | --version", "```", ""]
    with open(os.path.join(HERE, "API.md"), "w") as fh:
        fh.write("\n".join(out))
    print(f"wrote docs/API.md ({len(fm.METRICS)} metrics)")


if __name__ == "__main__":
    main()
