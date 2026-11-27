# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Human-perception-inspired metrics. Argument order (A, B, F).

Q_CV.m  Chen & Varshney 2007 (contrast-sensitivity-filtered difference energy; lower is better)
Q_CB.m  Chen & Blum 2009 (local contrast preservation with a CSF band-pass filter)
"""
import numpy as np

from ._matlab import as_double, blkproc, filter2_same, freqspace_meshgrid, im2double, mround


def _normalize1(x):
    """normalize1.m: min-max to 0..255 and MATLAB round (data unchanged when it is all zero)."""
    x = as_double(x)
    da, xiao = x.max(), x.min()
    if da == 0 and xiao == 0:
        return x
    return mround((x - xiao) / (da - xiao) * 255)


def _csf_filter(img, H):
    """ifft2(ifftshift(fftshift(fft2(img)) .* H)) (complex result, as in MATLAB)."""
    return np.fft.ifft2(np.fft.ifftshift(np.fft.fftshift(np.fft.fft2(img)) * H))


# ---------------------------------------------------------------- Q_CV.m (Chen & Varshney 2007)
def Q_CV(A, B, F):
    """Chen-Varshney metric (Q_CV.m): CSF-filtered squared differences between each source and F
    in 16x16 blocks, weighted by the block's Sobel edge energy (alpha = 5). Lower is better.

    In MATLAB the result is complex with an imaginary part at round-off level (~1e-14, the CSF is
    not Hermitian-symmetric on the even frequency grid); the real part is returned."""
    windowSize, alpha = 16, 5
    im1, im2, fim = _normalize1(A), _normalize1(B), _normalize1(F)
    flt1 = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], dtype=np.float64)
    flt2 = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]], dtype=np.float64)

    def grad(x):
        gx, gy = filter2_same(flt1, x), filter2_same(flt2, x)
        return np.sqrt(gx * gx + gy * gy)
    hang, lie = im1.shape
    ramda1 = blkproc(grad(im1), windowSize, lambda x: np.sum(x ** alpha))
    ramda2 = blkproc(grad(im2), windowSize, lambda x: np.sum(x ** alpha))
    u, v = freqspace_meshgrid(hang, lie)
    u, v = lie / 8 * u, hang / 8 * v
    r = np.sqrt(u ** 2 + v ** 2)
    theta_m = 2.6 * (0.0192 + 0.144 * r) * np.exp(-(0.144 * r) ** 1.1)
    Df1 = _csf_filter(im1 - fim, theta_m)
    Df2 = _csf_filter(im2 - fim, theta_m)
    D1 = blkproc(Df1, windowSize, lambda x: np.mean(x ** 2))
    D2 = blkproc(Df2, windowSize, lambda x: np.mean(x ** 2))
    Q = np.sum(ramda1 * D1 + ramda2 * D2) / np.sum(ramda1 + ramda2)
    return float(np.real(Q))


# ---------------------------------------------------------------- Q_CB.m (Chen & Blum 2009)
def _gaussian31(sigma):
    """gaussian2d in Q_CB.m: a 31x31 Gaussian (not normalised to unit sum), independent of the
    image size arguments."""
    x, y = np.meshgrid(np.arange(-15, 16), np.arange(-15, 16))
    return np.exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * np.pi * sigma * sigma)


def Q_CB(A, B, F):
    """Chen-Blum metric (Q_CB.m): preservation of CSF-filtered local contrast, k=h=1, p=3, q=2,
    Z=1e-4, Gaussians sigma 2 and 4. Higher is better."""
    f0, f1, a = 15.3870, 1.3456, 0.7622
    k, h, p, q, Z = 1, 1, 3, 2, 0.0001
    im1, im2, fim = (_normalize1(im2double(X)) for X in (A, B, F))
    hang, lie = im1.shape
    u, v = freqspace_meshgrid(hang, lie)
    u, v = (lie / 30) * u, (hang / 30) * v
    r = np.sqrt(u ** 2 + v ** 2)
    Sd = np.exp(-(r / f0) ** 2) - a * np.exp(-(r / f1) ** 2)
    G1, G2 = _gaussian31(2), _gaussian31(4)

    def cp(img):
        c = np.abs(filter2_same(G1, img) / filter2_same(G2, img) - 1)
        return (k * c ** p) / (h * c ** q + Z)
    C1P, C2P, CfP = cp(_csf_filter(im1, Sd)), cp(_csf_filter(im2, Sd)), cp(_csf_filter(fim, Sd))

    def pres(CP):
        mask = (CP < CfP).astype(np.float64)
        return (CP / CfP) * mask + (CfP / CP) * (1 - mask)
    Q1F, Q2F = pres(C1P), pres(C2P)
    den = C1P * C1P + C2P * C2P
    Q = (C1P * C1P) / den * Q1F + (C2P * C2P) / den * Q2F
    return float(np.mean(Q))


__all__ = ["Q_CV", "Q_CB"]
