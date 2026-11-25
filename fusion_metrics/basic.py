# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Line-by-line ports of Q_MLI.m, Q_AG.m, Q_En.m, Q_SF.m, Q_PSNR.m, Q_CC.m, Q_SCD.m, Q_MI.m.

Inputs follow the MATLAB files: grayscale images, uint8 or double with values in [0, 255]
(Q_MLI, Q_AG, Q_En also accept [0, 1]). Argument order is the MATLAB order: (A, B, F).
"""
import numpy as np

from ._matlab import EPS, as_double, mround


def _unit_range(F, name):
    F = as_double(F)
    if F.min() >= 0 and F.max() <= 1:
        return F
    if F.min() >= 0 and F.max() <= 255:
        return F / 255
    raise ValueError(f"{name}: input image must have values in either [0,1] or [0,255].")


def Q_MLI(F):
    """Mean luminance intensity on [0, 1] (Q_MLI.m)."""
    return float(np.mean(_unit_range(F, "Q_MLI")))


def Q_AG(F):
    """Average gradient on the [0, 1] scale (Q_AG.m: imgradientxy 'intermediate')."""
    I = _unit_range(F, "Q_AG")
    Gx = np.zeros_like(I)
    Gy = np.zeros_like(I)
    Gx[:, :-1] = I[:, 1:] - I[:, :-1]
    Gy[:-1, :] = I[1:, :] - I[:-1, :]
    return float(np.sum(np.sqrt(Gx * Gx + Gy * Gy)) / Gx.size)


def Q_En(F):
    """Shannon entropy in bits (Q_En.m: MATLAB entropy(), i.e. im2uint8 + 256-bin imhist)."""
    I = _unit_range(F, "Q_En")
    u = np.clip(np.floor(I * 255 + 0.5), 0, 255).astype(np.int64)   # im2uint8
    p = np.bincount(u.ravel(), minlength=256).astype(np.float64)
    p = p[p != 0] / u.size
    return float(-np.sum(p * np.log2(p)))


# ---------------------------------------------------------------- spatial frequency (Q_SF.m)
def _sf_components(im):
    im = as_double(im)
    hang, lie = im.shape
    RF = np.sqrt(np.sum((im[:, :-1] - im[:, 1:]) ** 2) / (hang * lie))
    CF = np.sqrt(np.sum((im[:-1, :] - im[1:, :]) ** 2) / (hang * lie))
    Wd = np.sqrt(2) / 2
    MDF = np.sqrt(Wd * np.sum((im[1:, 1:] - im[:-1, :-1]) ** 2) / (hang * lie))
    SDF = np.sqrt(Wd * np.sum((im[1:, :-1] - im[:-1, 1:]) ** 2) / (hang * lie))
    return RF, CF, MDF, SDF


def SF(F):
    """Four-directional spatial frequency of F: the SFf term inside Q_SF.m (Zheng et al. 2007)."""
    RF, CF, MDF, SDF = _sf_components(F)
    return float(np.sqrt(RF * RF + CF * CF + MDF * MDF + SDF * SDF))


def _sf_reference(im1, im2):
    a, b = as_double(im1), as_double(im2)
    hang, lie = a.shape
    Wd = np.sqrt(2) / 2

    def g(d1, d2):
        m = np.maximum(np.abs(d1), np.abs(d2))
        return np.sum(m * m)
    RFr = g(a[:, :-1] - a[:, 1:], b[:, :-1] - b[:, 1:]) / (hang * lie)
    CFr = g(a[:-1, :] - a[1:, :], b[:-1, :] - b[1:, :]) / (hang * lie)
    MDFr = Wd * g(a[1:, 1:] - a[:-1, :-1], b[1:, 1:] - b[:-1, :-1]) / (hang * lie)
    SDFr = Wd * g(a[1:, :-1] - a[:-1, 1:], b[1:, :-1] - b[:-1, 1:]) / (hang * lie)
    return float(np.sqrt(RFr + CFr + MDFr + SDFr))


def Q_SF(A, B, F):
    """Ratio of spatial-frequency error rSFe = (SF_F - SF_R) / SF_R (Q_SF.m)."""
    SFf = SF(F)
    SFr = _sf_reference(A, B)
    return 0.0 if SFr < 1e-10 else float((SFf - SFr) / SFr)


# ---------------------------------------------------------------- reference-based
def Q_PSNR(A, B, F):
    """Mean of psnr(F, A, 255) and psnr(F, B, 255) (Q_PSNR.m)."""
    f = as_double(F)

    def one(X):
        mse = np.mean((f - as_double(X)) ** 2)
        return np.inf if mse == 0 else 10 * np.log10(255.0 ** 2 / mse)
    return float((one(A) + one(B)) / 2)


def _corr(x, y):
    x = as_double(x).ravel(); y = as_double(y).ravel()
    xc = x - x.mean(); yc = y - y.mean()
    return np.sum(xc * yc) / np.sqrt(np.sum(xc * xc) * np.sum(yc * yc))


def Q_CC(A, B, F):
    """Mean correlation coefficient (r_AF + r_BF) / 2 (Q_CC.m, corrcoef)."""
    return float((_corr(A, F) + _corr(B, F)) / 2)


def _cc_eps(X, Y):
    X = as_double(X); Y = as_double(Y)
    Xc = X - X.mean(); Yc = Y - Y.mean()
    den = np.sqrt(np.sum(Xc * Xc) * np.sum(Yc * Yc))
    return 0.0 if den < EPS else float(np.sum(Xc * Yc) / den)


def Q_SCD(A, B, F):
    """Sum of the correlations of differences (Q_SCD.m): corr(F-B, A) + corr(F-A, B)."""
    A, B, F = as_double(A), as_double(B), as_double(F)
    return _cc_eps(F - B, A) + _cc_eps(F - A, B)


# ---------------------------------------------------------------- mutual information (Q_MI.m)
def _levels(X, uint8_input):
    """Histogram index used by Q_MI.m (grey_matrix + 1). With uint8 input MATLAB saturates
    255 + 1 to 255, so grey level 255 is counted in the bin of grey level 254."""
    X = np.asarray(X)
    if uint8_input:
        return np.minimum(X.astype(np.int64) + 1, 255) - 1
    Xd = X.astype(np.float64)
    if np.any(Xd != np.floor(Xd)) or Xd.min() < 0 or Xd.max() > 255:
        raise ValueError("Q_MI: grey levels must be integers in [0, 255]")
    return Xd.astype(np.int64)


def _H(counts):
    p = counts / counts.sum()
    p = p[p != 0]
    return float(np.sum(-p * np.log2(p)))


def Q_MI(A, B, F, grey_level=256, uint8_input=None):
    """MI(A,F) + MI(B,F) in bits (Q_MI.m).

    uint8_input reproduces MATLAB's saturating uint8 arithmetic in `grey_matrix + 1`:
    None (default) follows the dtype of each array, exactly as MATLAB would.
    """
    idx = {}
    for k, X in (("A", A), ("B", B), ("F", F)):
        u8 = (np.asarray(X).dtype == np.uint8) if uint8_input is None else uint8_input
        idx[k] = _levels(X, u8)
    HA = _H(np.bincount(idx["A"].ravel(), minlength=grey_level).astype(float))
    HB = _H(np.bincount(idx["B"].ravel(), minlength=grey_level).astype(float))
    HF = _H(np.bincount(idx["F"].ravel(), minlength=grey_level).astype(float))

    def Hab(X, Y):
        j = np.bincount((X.ravel() * 256 + Y.ravel()), minlength=256 * 256).astype(float)
        return _H(j)
    return float((HA + HF - Hab(idx["F"], idx["A"])) + (HB + HF - Hab(idx["F"], idx["B"])))


__all__ = ["Q_MLI", "Q_AG", "Q_En", "SF", "Q_SF", "Q_PSNR", "Q_CC", "Q_SCD", "Q_MI", "mround"]
