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
from .basic import Q_MLI, Q_AG, Q_En, SF, Q_SF, Q_PSNR, Q_CC, Q_SCD, Q_MI
from .info import Q_ABF, Q_VIFF, Q_M, Q_NCIE, Q_TE_Nava
from .fsim import Q_FSIM
from .fmi import Q_FMI, edge_sobel
from .variants import AG_VIFB, VIF_Sheikh

__version__ = "1.0.0"
__author__ = "Vu Viet Anh (Vũ Việt Anh)"
__email__ = "vietanhresearcher@gmail.com"

# name -> (function, needs source images)
METRICS = {
    "Q_MLI": (Q_MLI, False), "Q_AG": (Q_AG, False), "Q_En": (Q_En, False), "SF": (SF, False),
    "Q_SF": (Q_SF, True), "Q_PSNR": (Q_PSNR, True), "Q_CC": (Q_CC, True), "Q_SCD": (Q_SCD, True),
    "Q_MI": (Q_MI, True), "Q_ABF": (Q_ABF, True), "Q_VIFF": (Q_VIFF, True), "Q_M": (Q_M, True),
    "Q_NCIE": (Q_NCIE, True), "Q_TE_Nava": (Q_TE_Nava, True), "Q_FSIM": (Q_FSIM, True),
    "Q_FMI": (Q_FMI, True),
    # widely used alternative definitions; MATLAB references in ../matlab_added/
    "AG_VIFB": (AG_VIFB, False), "VIF_Sheikh": (VIF_Sheikh, True),
}

# Ready-made metric groups (Liu et al., TPAMI 2012 taxonomy, plus image statistics).
PRESETS = {
    "image": ["Q_MLI", "Q_En", "SF", "AG_VIFB"],                       # statistics of F alone
    "information": ["Q_MI", "Q_NCIE", "Q_TE_Nava", "Q_FMI"],           # information theory
    "feature": ["Q_ABF", "Q_M", "Q_SF", "Q_FSIM"],                     # edges / features
    "similarity": ["Q_PSNR", "Q_CC", "Q_SCD", "VIF_Sheikh", "Q_VIFF"],  # similarity to sources
}
PRESETS["all"] = list(METRICS)


def evaluate(A, B, F, names=None):
    """Compute the requested metrics (a list of names or a PRESETS key; default: all)."""
    if isinstance(names, str):
        names = PRESETS[names]
    names = names or list(METRICS)
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

__all__ = list(METRICS) + ["METRICS", "PRESETS", "evaluate", "evaluate_files", "read_gray",
                           "read_gray_cv", "rgb2gray", "edge_sobel", "__version__"]
