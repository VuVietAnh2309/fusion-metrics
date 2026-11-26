# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Ports of Q_ABF.m, Q_VIFF.m, Q_M.m, Q_NCIE.m, Q_NMI.m, Q_CE.m, Q_TE_Nava.m. Argument order (A, B, F)."""
import numpy as np

from ._matlab import as_double, conv2_same, filter2_valid, fspecial_gaussian, mround, to_uint8


# ---------------------------------------------------------------- Q_ABF.m (Petrovic Q^{AB/F})
def _sobel_ga(p):
    h1 = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]], dtype=np.float64)
    h3 = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], dtype=np.float64)
    Sx = conv2_same(p, h3)
    Sy = conv2_same(p, h1)
    g = np.sqrt(Sx ** 2 + Sy ** 2)
    a = np.full_like(Sx, np.pi / 2)
    nz = Sx != 0
    a[nz] = np.arctan(Sy[nz] / Sx[nz])
    return g, a


def Q_ABF(A, B, F):
    """Gradient-based fusion quality Q^{AB/F} (Q_ABF.m)."""
    Tg, kg, Dg = 0.9994, -15, 0.5
    Ta, ka, Da = 0.9879, -22, 0.8
    gA, aA = _sobel_ga(as_double(A))
    gB, aB = _sobel_ga(as_double(B))
    gF, aF = _sobel_ga(as_double(F))

    def q(gS, aS):
        G = np.empty_like(gS)
        m1 = gS > gF
        m2 = gS == gF
        m3 = ~(m1 | m2)
        G[m1] = gF[m1] / gS[m1]
        G[m2] = gF[m2]
        G[m3] = gS[m3] / gF[m3]
        Ag = 1 - np.abs(aS - aF) / (np.pi / 2)
        return (Tg / (1 + np.exp(kg * (G - Dg)))) * (Ta / (1 + np.exp(ka * (Ag - Da))))
    QAF, QBF = q(gA, aA), q(gB, aB)
    deno = np.sum(gA + gB)
    return 0.0 if deno == 0 else float(np.sum(QAF * gA + QBF * gB) / deno)


# ---------------------------------------------------------------- Q_VIFF.m (Han et al. 2013)
def _vid_vind_g(ref, dist, sq):
    num, den, G = [], [], []
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
        G.append(g)
        num.append(np.log10(1 + g ** 2 * s1 / (sv + sq)))
        den.append(np.log10(1 + s1 / sq))
    return num, den, G


def Q_VIFF(A, B, F):
    """Visual information fidelity for fusion, VIFF (Han et al., 2013; Q_VIFF.m), grayscale."""
    sq = 0.005 * 255 * 255
    C = 1e-7
    T1N, T1D, T1G = _vid_vind_g(A, F, sq)
    T2N, T2D, T2G = _vid_vind_g(B, F, sq)
    p = np.array([1, 0, 0.15, 1]) / 2.15
    out = np.zeros(4)
    for i in range(4):
        L = T1G[i] < T2G[i]
        Z12 = np.where(L, T1N[i], T2N[i])
        M12 = np.where(L, T1D[i], T2D[i])
        VID = np.sum(Z12 + C)
        VIND = np.sum(M12 + C)
        out[i] = 0 if VIND < 1e-10 else VID / VIND
    return float(np.sum(out * p))


# ---------------------------------------------------------------- Q_M.m (Wang & Liu, Haar)
def _dwt2_haar(x):
    """MATLAB dwt2(x, 'haar') with the default 'sym' extension -> (cA, cH, cV, cD)."""
    x = as_double(x)

    def ext(v, axis):   # half-point symmetric extension by one sample when the length is odd
        if v.shape[axis] % 2:
            last = np.take(v, [-1], axis=axis)
            v = np.concatenate([v, last], axis=axis)
        return v
    x = ext(ext(x, 0), 1)
    s = 1 / np.sqrt(2)
    lo_c = (x[:, 0::2] + x[:, 1::2]) * s          # columns low-pass (along rows)
    hi_c = (x[:, 0::2] - x[:, 1::2]) * s
    cA = (lo_c[0::2] + lo_c[1::2]) * s
    cH = (lo_c[0::2] - lo_c[1::2]) * s
    cV = (hi_c[0::2] + hi_c[1::2]) * s
    cD = (hi_c[0::2] - hi_c[1::2]) * s
    return cA, cH, cV, cD


def Q_M(A, B, F):
    """Multiscale wavelet edge-preservation metric Q_M (Q_M.m, 2-level Haar)."""
    alpha = np.array([2 / 3, 1 / 3])
    L1, L2, Lf = as_double(A), as_double(B), as_double(F)
    Q = np.zeros(2)
    for s in range(2):
        L1, H1, V1, D1 = _dwt2_haar(L1)
        L2, H2, V2, D2 = _dwt2_haar(L2)
        Lf, Hf, Vf, Df = _dwt2_haar(Lf)
        EPA = (np.exp(-np.abs(H1 - Hf)) + np.exp(-np.abs(V1 - Vf)) + np.exp(-np.abs(D1 - Df))) / 3
        EPB = (np.exp(-np.abs(H2 - Hf)) + np.exp(-np.abs(V2 - Vf)) + np.exp(-np.abs(D2 - Df))) / 3
        wA = H1 ** 2 + V1 ** 2 + D1 ** 2
        wB = H2 ** 2 + V2 ** 2 + D2 ** 2
        Q[s] = np.sum(EPA * wA + EPB * wB) / np.sum(wA + wB)
    return float(np.prod(Q ** alpha))


# ---------------------------------------------------------------- Q_NCIE.m
def _normalize1(x):
    x = as_double(x)
    da, xiao = x.max(), x.min()
    if da == 0 and xiao == 0:
        return x
    return mround((x - xiao) / (da - xiao) * 255)


def _ncc(a, b):
    b_ = 256
    h = np.bincount((a.astype(np.int64).ravel() * 256 + b.astype(np.int64).ravel()),
                    minlength=256 * 256).reshape(256, 256).astype(np.float64)
    h = h / h.sum()
    m1 = h.sum(axis=0)
    m2 = h.sum(axis=1)
    Hx = -np.sum(m1 * np.log2(m1 + (m1 == 0)))
    Hy = -np.sum(m2 * np.log2(m2 + (m2 == 0)))
    Hxy = -np.sum(h * np.log2(h + (h == 0)))
    return Hx / np.log2(b_) + Hy / np.log2(b_) - Hxy / np.log2(b_)


def Q_NCIE(A, B, F):
    """Nonlinear correlation information entropy (Q_NCIE.m)."""
    a, b, f = _normalize1(A), _normalize1(B), _normalize1(F)
    K = 3
    xy, xf, yf = _ncc(a, b), _ncc(a, f), _ncc(b, f)
    R = np.array([[1, xy, xf], [xy, 1, yf], [xf, yf, 1]])
    r = np.linalg.eigvalsh(R)
    HR = -np.sum(r * np.log2(r / K) / K) / np.log2(256)
    return float(1 - HR)


# ---------------------------------------------------------------- Q_NMI.m (Hossny et al.)
def _mutual_info(im1, im2):
    """mutual_info.m (Z. Liu): joint histogram h(im1+1, im2+1) -> (MI, H_xy, H_x, H_y).
    As in the .m file, H_x is computed from the column sums (the marginal of im2) and H_y from
    the row sums (the marginal of im1); Q_NMI only uses their sum, so the naming does not matter."""
    h = np.bincount((im1.astype(np.int64).ravel() * 256 + im2.astype(np.int64).ravel()),
                    minlength=256 * 256).reshape(256, 256).astype(np.float64)
    h = h / h.sum()
    m1 = h.sum(axis=0)
    m2 = h.sum(axis=1)
    H_x = -np.sum(m1 * np.log2(m1 + (m1 == 0)))
    H_y = -np.sum(m2 * np.log2(m2 + (m2 == 0)))
    H_xy = -np.sum(h * np.log2(h + (h == 0)))
    return H_x + H_y - H_xy, H_xy, H_x, H_y


def Q_NMI(A, B, F):
    """Normalized mutual information 2*[MI(A,F)/(H(A)+H(F)) + MI(B,F)/(H(B)+H(F))]
    (Q_NMI.m; inputs min-max normalised and rounded to 0..255 by normalize1)."""
    a, b, f = _normalize1(A), _normalize1(B), _normalize1(F)
    I_fx, _, H_x, H_f1 = _mutual_info(a, f)
    I_fy, _, H_y, H_f2 = _mutual_info(b, f)
    return float(2 * (I_fx / (H_f1 + H_x) + I_fy / (H_f2 + H_y)))


# ---------------------------------------------------------------- Q_CE.m (cross entropy)
def _prepare_image(img):
    """prepare_image in Q_CE.m: uint8 is kept; otherwise double, x255 when max <= 1, clipped to
    [0, 255] and converted with uint8() (round halves away from zero)."""
    img = np.asarray(img)
    if img.dtype == np.uint8:
        return img
    img = img.astype(np.float64)
    if img.max() <= 1:
        img = img * 255
    return to_uint8(np.clip(img, 0, 255))


def _kl(hP, hF):
    mask = hP > 0
    return float(np.sum(hP[mask] * np.log2(hP[mask] / np.maximum(hF[mask], 1e-12))))


def Q_CE(A, B, F):
    """Cross entropy 0.5*[D_KL(h_A||h_F) + D_KL(h_B||h_F)] of the 256-bin histograms
    (Q_CE.m; lower is better)."""
    hA, hB, hF = (np.bincount(_prepare_image(X).ravel(), minlength=256).astype(np.float64)
                  for X in (A, B, F))
    hA, hB, hF = hA / hA.sum(), hB / hB.sum(), hF / hF.sum()
    return 0.5 * (_kl(hA, hF) + _kl(hB, hF))


# ---------------------------------------------------------------- Q_TE_Nava.m (Tsallis)
def _tsallis_entropy(X, q):
    c = np.bincount(X.ravel(), minlength=256).astype(np.float64)
    p = c / c.sum()
    p = p[p > 0]
    return (1 / (1 - q)) * (np.sum(p ** q) - 1)


def _tsallis_mi(X, Y, q):
    j = np.bincount(X.ravel() * 256 + Y.ravel(), minlength=256 * 256).reshape(256, 256).astype(np.float64)
    pxy = j / j.sum()
    px = pxy.sum(axis=1)
    py = pxy.sum(axis=0)
    i, k = np.nonzero(pxy.T)          # MATLAB find() is column-major; order only affects rounding
    k, i = i, k
    p = pxy[i, k]
    return (1 / (1 - q)) * (1 - np.sum((p ** q) / ((px[i] * py[k]) ** (q - 1))))


def Q_TE_Nava(A, B, F, q=0.43137):
    """Tsallis-entropy fusion metric, Nava et al. normalisation (Q_TE_Nava.m)."""
    a, b, f = (to_uint8(as_double(X)).astype(np.int64) for X in (A, B, F))
    I_FA = _tsallis_mi(f, a, q)
    I_FB = _tsallis_mi(f, b, q)
    H_A = _tsallis_entropy(a, q)
    H_B = _tsallis_entropy(b, q)
    I_AB = _tsallis_mi(a, b, q)
    return float((I_FA + I_FB) / (H_A + H_B - I_AB))


__all__ = ["Q_ABF", "Q_VIFF", "Q_M", "Q_NCIE", "Q_NMI", "Q_CE", "Q_TE_Nava"]
