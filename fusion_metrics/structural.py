# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Ports of the SSIM-family metrics. Argument order (A, B, F).

MS_SSIM.m  SSIM (Z. Wang 2003), mean over the two sources
Q_UIQI.m   UIQI (Z. Wang & A. C. Bovik 2001), mean over the two sources
Q_S.m      Piella & Heijmans 2003: Q_S (sw=1), weighted Q_W (sw=2), edge-dependent Q_E (sw=3)
Q_C.m      Cvejic et al. 2005 (covariance-weighted SSIM)
Q_Y.m      Yang et al. 2008
Q_H.m      local-entropy weighted SSIM of Hossny-type (7x7 windows)
"""
import numpy as np

from ._matlab import as_double, filter2_same, filter2_valid, fspecial_gaussian, mround, mtimes_elem


# ---------------------------------------------------------------- ssim_index.m / ssim_yang
def _flat_windows(img, w):
    """True for every 'valid' w x w window in which img is constant."""
    from numpy.lib.stride_tricks import sliding_window_view
    v = sliding_window_view(img, (w, w))
    return v.max(axis=(-1, -2)) == v.min(axis=(-1, -2))


def _ssim_parts(img1, img2, wsize, sigma, C1, C2, exact_flat=False):
    """Shared body of ssim_index.m and ssim_yang: Gaussian window, 'valid' filtering.
    Returns (ssim_map, sigma1_sq, sigma2_sq, sigma12)."""
    M, N = np.shape(img1)
    if M < 11 or N < 11:
        # the .m files return without assigning their outputs, which is an error in MATLAB
        raise ValueError("SSIM-based metrics need images of at least 11x11")
    window = fspecial_gaussian(wsize, sigma)
    window = window / window.sum()
    img1, img2 = as_double(img1), as_double(img2)
    mu1 = filter2_valid(window, img1)
    mu2 = filter2_valid(window, img2)
    mu1_sq, mu2_sq, mu1_mu2 = mu1 * mu1, mu2 * mu2, mu1 * mu2
    sigma1_sq = filter2_valid(window, img1 * img1) - mu1_sq
    sigma2_sq = filter2_valid(window, img2 * img2) - mu2_sq
    sigma12 = filter2_valid(window, img1 * img2) - mu1_mu2
    if exact_flat:
        # Where an image is constant over the window its local variance (and the covariance) is
        # exactly 0, but filter2 leaves round-off of order 1e-11. With C ~ 2e-16 (ssim_yang) the .m
        # result then depends on the summation order (MATLAB, Octave and NumPy all differ, e.g.
        # 0.67 / 1.5e5 for a window whose exact value is 0.9986). These windows are evaluated in
        # exact arithmetic; everywhere else the .m computation is unchanged.
        f1, f2 = _flat_windows(img1, wsize), _flat_windows(img2, wsize)
        sigma1_sq = np.where(f1, 0.0, sigma1_sq)
        sigma2_sq = np.where(f2, 0.0, sigma2_sq)
        sigma12 = np.where(f1 | f2, 0.0, sigma12)
    ssim_map = ((2 * mu1_mu2 + C1) * (2 * sigma12 + C2)) / ((mu1_sq + mu2_sq + C1) * (sigma1_sq + sigma2_sq + C2))
    return ssim_map, sigma1_sq, sigma2_sq, sigma12


def _ssim_index_map(img1, img2):
    """ssim_index.m with its default arguments (11x11, sigma 1.5, K=[0.01 0.03], L=255)."""
    return _ssim_parts(img1, img2, 11, 1.5, (0.01 * 255) ** 2, (0.03 * 255) ** 2)


def _ssim_yang(img1, img2):
    """ssim_yang in Q_Y.m / Q_C.m: 7x7 Gaussian window (sigma 1.5), C1 = C2 = 2e-16."""
    return _ssim_parts(img1, img2, 7, 1.5, 2e-16, 2e-16, exact_flat=True)


# ---------------------------------------------------------------- MS_SSIM.m
def _ssim_index(img1, img2):
    return float(np.mean(_ssim_index_map(img1, img2)[0]))


def MS_SSIM(A, B, F):
    """Mean SSIM (ssim_index(A,F) + ssim_index(B,F)) / 2 (MS_SSIM.m). Despite its name this is
    single-scale SSIM averaged over the two sources, exactly as in the reference file."""
    return float((_ssim_index(A, F) + _ssim_index(B, F)) / 2)


# ---------------------------------------------------------------- Q_UIQI.m (UIQI)
def _uiqi(img1, img2, block_size=8):
    """UIQI.m (img_qi) with the default 8x8 sliding window; returns mean2(quality_map).

    Like the .m file, the inputs are NOT converted to double: for uint8 images the products
    img1.*img1, img2.*img2 and img1.*img2 are uint8 in MATLAB, i.e. rounded and saturated at 255
    (mtimes_elem reproduces this). filter2 then works in double."""
    img1, img2 = np.asarray(img1), np.asarray(img2)
    N = block_size ** 2
    k = np.ones((block_size, block_size))
    img1_sq = mtimes_elem(img1, img1)
    img2_sq = mtimes_elem(img2, img2)
    img12 = mtimes_elem(img1, img2)
    img1_sum = filter2_valid(k, as_double(img1))
    img2_sum = filter2_valid(k, as_double(img2))
    img1_sq_sum = filter2_valid(k, as_double(img1_sq))
    img2_sq_sum = filter2_valid(k, as_double(img2_sq))
    img12_sum = filter2_valid(k, as_double(img12))
    img12_sum_mul = img1_sum * img2_sum
    img12_sq_sum_mul = img1_sum * img1_sum + img2_sum * img2_sum
    numerator = 4 * (N * img12_sum - img12_sum_mul) * img12_sum_mul
    denominator1 = N * (img1_sq_sum + img2_sq_sum) - img12_sq_sum_mul
    denominator = denominator1 * img12_sq_sum_mul
    quality_map = np.ones(denominator.shape)
    idx = (denominator1 == 0) & (img12_sq_sum_mul != 0)
    quality_map[idx] = 2 * img12_sum_mul[idx] / img12_sq_sum_mul[idx]
    idx = denominator != 0
    quality_map[idx] = numerator[idx] / denominator[idx]
    return float(np.mean(quality_map))


def Q_UIQI(A, B, F):
    """Mean universal image quality index (UIQI(A,F) + UIQI(B,F)) / 2 (Q_UIQI.m).

    Pass double images to get the textbook UIQI; with uint8 inputs the result follows MATLAB's
    uint8 arithmetic inside UIQI.m (saturated squares), exactly as the reference file does."""
    return float((_uiqi(A, F) + _uiqi(B, F)) / 2)


# ---------------------------------------------------------------- Q_S.m (Piella & Heijmans 2003)
def _two_plane_max_weights(s1, s2):
    """[Cw, U] = max(cat(3, s1, s2), [], 3); cw = Cw / sum(sum(Cw))."""
    Cw = np.maximum(s1, s2)
    return Cw / Cw.sum()


def _lambda(s1, s2):
    """The 0.5 guard of Q_S.m / Q_Y.m: where s1 + s2 == 0 both get +0.5, then s1 / (s1 + s2)."""
    test = (s1 + s2 == 0) * 0.5
    s1, s2 = s1 + test, s2 + test
    return s1 / (s1 + s2), s1, s2


def Q_S(A, B, F, sw=1):
    """Piella's metric (Q_S.m = metricPeilla(img1, img2, fuse, sw)).
    sw=1: Q_S, sw=2: weighted Q_W, sw=3: edge-dependent Q_E (alpha = 1)."""
    img1, img2, fuse = as_double(A), as_double(B), as_double(F)
    if sw == 3:
        flt1 = np.array([[1, 0, -1], [1, 0, -1], [1, 0, -1]], dtype=np.float64)
        flt2 = np.array([[1, 1, 1], [0, 0, 0], [-1, -1, -1]], dtype=np.float64)

        def grad(x):
            gx, gy = filter2_same(flt1, x), filter2_same(flt2, x)
            return np.sqrt(gx * gx + gy * gy)
        img1, img2, fuse = grad(img1), grad(img2), grad(fuse)
    _, s1, s2, _ = _ssim_index_map(img1, img2)
    ramda, s1, s2 = _lambda(s1, s2)
    m1 = _ssim_index_map(fuse, img1)[0]
    m2 = _ssim_index_map(fuse, img2)[0]
    Q = ramda * m1 + (1 - ramda) * m2
    if sw == 1:
        return float(np.mean(Q))
    if sw not in (2, 3):
        raise ValueError("sw must be 1, 2 or 3")
    Qw = float(np.sum(_two_plane_max_weights(s1, s2) * Q))
    return Qw if sw == 2 else Qw * Qw ** 1          # Qe = Qw * Qw^alpha, alpha = 1


def Q_Sw(A, B, F):
    """Weighted fusion quality index Q_W of Piella (Q_S.m with sw = 2)."""
    return Q_S(A, B, F, sw=2)


def Q_Se(A, B, F):
    """Edge-dependent fusion quality index Q_E of Piella (Q_S.m with sw = 3, alpha = 1)."""
    return Q_S(A, B, F, sw=3)


# ---------------------------------------------------------------- Q_Y.m (Yang et al. 2008)
def Q_Y(A, B, F):
    """Yang et al.'s SSIM-based fusion metric (Q_Y.m)."""
    im1, im2, fim = as_double(A), as_double(B), as_double(F)
    map1, s1, s2, _ = _ssim_yang(im1, im2)
    map2 = _ssim_yang(im1, fim)[0]
    map3 = _ssim_yang(im2, fim)[0]
    bin_map = map1 >= 0.75
    ramda, _, _ = _lambda(s1, s2)
    Q1 = (ramda * map2 + (1 - ramda) * map3) * bin_map
    Q2 = np.maximum(map2, map3) * (~bin_map)
    return float(np.mean(Q1 + Q2))


# ---------------------------------------------------------------- Q_C.m (Cvejic et al. 2005)
def Q_C(A, B, F):
    """Cvejic et al.'s metric (Q_C.m): SSIM maps weighted by the local covariances sigma_XF,
    sigma_YF, with the weighted value forced to 0 / 1 where the weight leaves [0, 1]."""
    im1, im2, fim = as_double(A), as_double(B), as_double(F)
    map2, _, _, sXF = _ssim_yang(im1, fim)
    map3, _, _, sYF = _ssim_yang(im2, fim)
    sim_w, _, _ = _lambda(sXF, sYF)
    sim = sim_w * map2 + (1 - sim_w) * map3
    sim[sim_w < 0] = 0
    sim[sim_w > 1] = 1
    return float(np.mean(sim))


# ---------------------------------------------------------------- Q_H.m
def _local_entropy_u8(img, w=7):
    """local_entropy_u8 in Q_H.m: entropy (bits) of the uint8 levels in every w x w window
    ('valid' positions), with MATLAB's summation order over the occupied bins."""
    lvl = np.clip(mround(img), 0, 255).astype(np.int64)
    from numpy.lib.stride_tricks import sliding_window_view
    win = sliding_window_view(lvl, (w, w)).reshape(-1, w * w)
    n = win.shape[0]
    counts = np.zeros((n, 256))
    np.add.at(counts, (np.repeat(np.arange(n), w * w), win.ravel()), 1)
    p = counts / (w * w)
    terms = np.zeros_like(p)
    nz = p > 0
    terms[nz] = p[nz] * np.log2(p[nz])
    E = -np.cumsum(terms, axis=1)[:, -1]
    return E.reshape(lvl.shape[0] - w + 1, lvl.shape[1] - w + 1)


def Q_H(A, B, F):
    """Local-entropy weighted SSIM metric of Q_H.m (7x7 windows; returns mean(1 - E))."""
    w, half = 7, 3
    E1, E2, Ef = (_local_entropy_u8(as_double(X), w) for X in (A, B, F))
    p1, p2 = E1 / E1.max(), E2 / E2.max()
    window = fspecial_gaussian(7, 1.5)
    window = window / window.sum()

    def lssim(a, b):
        mu1, mu2 = filter2_valid(window, a), filter2_valid(window, b)
        s1 = filter2_valid(window, a * a) - mu1 * mu1
        s2 = filter2_valid(window, b * b) - mu2 * mu2
        s12 = filter2_valid(window, a * b) - mu1 * mu2
        return ((2 * mu1 * mu2 + 2e-16) * (2 * s12 + 2e-16)) / ((mu1 * mu1 + mu2 * mu2 + 2e-16) * (s1 + s2 + 2e-16))
    m1, m2 = lssim(E1, Ef), lssim(E2, Ef)

    def crop(x):
        return x[half:x.shape[0] - half, half:x.shape[1] - half]
    p1, p2, E1c, E2c, Efc = crop(p1), crop(p2), crop(E1), crop(E2), crop(Ef)
    S_xf = E1c / (Efc + 1e-12)
    S_yf = E2c / (Efc + 1e-12)
    ramda = S_xf / (S_xf + S_yf + 1e-12)
    E = ramda * (p1 * m1) + (1 - ramda) * (p2 * m2)
    return float(np.mean(1 - E))


__all__ = ["MS_SSIM", "Q_UIQI", "Q_S", "Q_Sw", "Q_Se", "Q_Y", "Q_C", "Q_H"]
