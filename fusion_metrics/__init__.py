# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""fusion_metrics: objective quality metrics for image fusion, faithfully ported from MATLAB.

Every function keeps the name and argument order of its reference .m file (../matlab/*.m):
Q_X(A, B, F) with A, B the two source images and F the fused image (grayscale, uint8 or double
in [0, 255]). No-reference metrics take only F.

Read images with `read_gray` (= MATLAB imread + rgb2gray) to get exactly the numbers the
MATLAB reference implementations produce.
"""
from ._matlab import read_gray, read_gray_cv, rgb2gray
from .basic import Q_MLI, Q_AG, Q_En, Q_CI, Q_EI, SF, Q_SF, Q_PSNR, Q_CC, Q_SCD, Q_MI
from .info import Q_ABF, Q_NABF, Q_VIFF, Q_M, Q_NCIE, Q_NMI, Q_CE, Q_TE, Q_TE_Nava
from .structural import MS_SSIM, Q_UIQI, Q_S, Q_Sw, Q_Se, Q_Y, Q_C, Q_H
from .perceptual import Q_CV, Q_CB
from .phase import Q_P
from .fsim import Q_FSIM
from .fmi import Q_FMI, edge_sobel
from .variants import AG_VIFB, VIF_Sheikh

import warnings as _warnings

import numpy as _np

__version__ = "1.2.0"
__author__ = "Vu Viet Anh (Vũ Việt Anh)"
__email__ = "vietanhresearcher@gmail.com"

# name -> (function, needs source images)
METRICS = {
    "Q_MLI": (Q_MLI, False), "Q_AG": (Q_AG, False), "Q_En": (Q_En, False), "SF": (SF, False),
    "Q_SF": (Q_SF, True), "Q_PSNR": (Q_PSNR, True), "Q_CC": (Q_CC, True), "Q_SCD": (Q_SCD, True),
    "Q_MI": (Q_MI, True), "Q_ABF": (Q_ABF, True), "Q_VIFF": (Q_VIFF, True), "Q_M": (Q_M, True),
    "Q_NCIE": (Q_NCIE, True), "Q_TE_Nava": (Q_TE_Nava, True), "Q_FSIM": (Q_FSIM, True),
    "Q_FMI": (Q_FMI, True),
    "Q_CI": (Q_CI, False), "Q_EI": (Q_EI, False), "Q_CE": (Q_CE, True), "Q_NMI": (Q_NMI, True),
    "MS_SSIM": (MS_SSIM, True), "Q_UIQI": (Q_UIQI, True),
    "Q_TE": (Q_TE, True), "Q_NABF": (Q_NABF, True), "Q_P": (Q_P, True),
    "Q_S": (Q_S, True), "Q_Sw": (Q_Sw, True), "Q_Se": (Q_Se, True), "Q_C": (Q_C, True),
    "Q_Y": (Q_Y, True), "Q_H": (Q_H, True), "Q_CV": (Q_CV, True), "Q_CB": (Q_CB, True),
    # widely used alternative definitions; MATLAB references in ../matlab_added/
    "AG_VIFB": (AG_VIFB, False), "VIF_Sheikh": (VIF_Sheikh, True),
}

# Ready-made metric groups (Liu et al., TPAMI 2012 taxonomy, plus image statistics).
PRESETS = {
    "image": ["Q_MLI", "Q_En", "Q_CI", "SF", "AG_VIFB", "Q_EI"],       # statistics of F alone
    "information": ["Q_MI", "Q_NMI", "Q_NCIE", "Q_TE", "Q_TE_Nava", "Q_FMI", "Q_CE"],   # information theory
    "feature": ["Q_ABF", "Q_NABF", "Q_M", "Q_SF", "Q_P", "Q_FSIM"],    # edges / features
    "structural": ["MS_SSIM", "Q_UIQI", "Q_S", "Q_Sw", "Q_Se", "Q_C", "Q_Y"],   # structural similarity
    "perception": ["Q_CV", "Q_CB"],                                    # human-perception inspired
    # Liu et al., TPAMI 2012 (the twelve metrics of its comparison)
    "liu2012": ["Q_NMI", "Q_TE", "Q_NCIE", "Q_ABF", "Q_M", "Q_SF", "Q_P", "Q_S", "Q_C", "Q_Y", "Q_CV", "Q_CB"],
    "similarity": ["Q_PSNR", "Q_CC", "Q_SCD", "VIF_Sheikh", "Q_VIFF"],  # similarity to sources
    # the table most infrared-visible fusion papers report (CDDFuse, DDFM, ...):
    # EN, SD, SF, MI, SCD, VIF, Qabf, SSIM
    "ivif": ["Q_En", "Q_CI", "SF", "Q_MI", "Q_SCD", "VIF_Sheikh", "Q_ABF", "MS_SSIM"],
}
PRESETS["all"] = list(METRICS)


class ScaleMismatchWarning(UserWarning):
    """The images do not look like they are in the [0, 255] range the metrics expect."""


def check_range(*images):
    """Warn (ScaleMismatchWarning) when images look normalised to [0, 1].

    The reference metrics expect grey levels in [0, 255]. An image whose maximum is <= 1 is
    almost certainly a normalised float image; mixing it with [0, 255] images makes gradient
    and SSIM-type metrics collapse silently (Q_ABF 0.40 -> 0.0004). Returns True if all
    images look fine."""
    peaks = [float(_np.max(x)) for x in images]
    low = [i for i, m in enumerate(peaks) if m <= 1.0]
    if not low:
        return True
    labels = "ABF" if len(images) == 3 else [f"#{i}" for i in range(len(images))]
    which = ", ".join(labels[i] for i in low)
    if len(low) < len(images):
        msg = (f"image(s) {which} have max <= 1 while the others reach "
               f"{max(peaks):.0f}: mixed value ranges make most metrics meaningless. "
               "Scale every image to [0, 255] (e.g. x * 255) before evaluating.")
    else:
        msg = ("all images have max <= 1; the metrics expect grey levels in [0, 255]. "
               "Multiply by 255 to get values comparable with published tables.")
    _warnings.warn(msg, ScaleMismatchWarning, stacklevel=3)
    return False


def evaluate(A, B, F, names=None, check=True):
    """Compute the requested metrics (a list of names or a PRESETS key; default: all).

    With check=True (default) a ScaleMismatchWarning is emitted when an image looks
    normalised to [0, 1]; see check_range()."""
    if isinstance(names, str):
        names = PRESETS[names]
    names = names or list(METRICS)
    if check:
        check_range(A, B, F)
    out = {}
    for n in names:
        fn, src = METRICS[n]
        out[n] = fn(A, B, F) if src else fn(F)
    return out


def evaluate_files(path_a, path_b, path_f, names=None, reader=None):
    """Same as evaluate() but reading the three files.

    reader=read_gray (default) follows MATLAB (imread + rgb2gray); reader=read_gray_cv uses the
    OpenCV gray conversion common in Python fusion code bases."""
    reader = reader or read_gray
    return evaluate(reader(path_a), reader(path_b), reader(path_f), names)


read_gray_paper = read_gray_cv  # backwards-compatible alias

__all__ = list(METRICS) + ["METRICS", "PRESETS", "evaluate", "evaluate_files", "check_range",
                           "ScaleMismatchWarning", "read_gray", "read_gray_cv", "rgb2gray",
                           "edge_sobel", "__version__"]
