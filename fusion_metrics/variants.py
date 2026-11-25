# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Two widely reported alternative definitions that differ from the Q_AG / Q_VIFF reference
files. Their MATLAB references are in ../matlab_added/.

AG_VIFB(F)        = AG_VIFB.m    (VIFB metricsAvg_gradient.m), 0..255 scale
VIF_Sheikh(A,B,F) = VIF_Sheikh.m (vifp_mscale(A,F) + vifp_mscale(B,F), Sheikh & Bovik 2006)
"""
import numpy as np

from ._matlab import as_double, filter2_valid, fspecial_gaussian


def AG_VIFB(F):
    F = as_double(F)
    r, c = F.shape
    dzdy, dzdx = np.gradient(F)            # MATLAB gradient(): central / one-sided differences
    s = np.sqrt((dzdx ** 2 + dzdy ** 2) / 2)
    return float(np.sum(s) / ((r - 1) * (c - 1)))


def _vifp_mscale(ref, dist):
    sigma_nsq = 2.0
    num = den = 0.0
    ref = as_double(ref); dist = as_double(dist)
    for scale in range(1, 5):
        N = 2 ** (4 - scale + 1) + 1
        win = fspecial_gaussian(N, N / 5)
        if scale > 1:
            ref = filter2_valid(win, ref)[::2, ::2]
            dist = filter2_valid(win, dist)[::2, ::2]
        mu1 = filter2_valid(win, ref)
        mu2 = filter2_valid(win, dist)
        s1 = filter2_valid(win, ref * ref) - mu1 * mu1
        s2 = filter2_valid(win, dist * dist) - mu2 * mu2
        s12 = filter2_valid(win, ref * dist) - mu1 * mu2
        s1[s1 < 0] = 0
        s2[s2 < 0] = 0
        g = s12 / (s1 + 1e-10)
        sv = s2 - g * s12
        m = s1 < 1e-10
        g[m] = 0; sv[m] = s2[m]; s1[m] = 0
        m = s2 < 1e-10
        g[m] = 0; sv[m] = 0
        m = g < 0
        sv[m] = s2[m]; g[m] = 0
        sv[sv <= 1e-10] = 1e-10
        num += np.sum(np.log10(1 + g ** 2 * s1 / (sv + sigma_nsq)))
        den += np.sum(np.log10(1 + s1 / sigma_nsq))
    return num / den


def VIF_Sheikh(A, B, F):
    return float(_vifp_mscale(A, F) + _vifp_mscale(B, F))


__all__ = ["AG_VIFB", "VIF_Sheikh"]
