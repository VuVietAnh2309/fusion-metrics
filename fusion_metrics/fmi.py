# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Port of Q_FMI.m (feature mutual information, Haghighat et al.) with its default
feature='edge', w=3, including MATLAB's default edge() (Sobel, automatic threshold, thinning).

The .m file is reproduced literally, including one quirk: inside each window the fused
window fSub is normalised in the A-branch and that normalised window is what the B-branch
compares with bSub and normalises again.
"""
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from ._matlab import as_double


# ---------------------------------------------------------------- MATLAB edge(I) (Sobel)
def _imfilter_replicate(a, k):
    """MATLAB imfilter(a, k, 'replicate') for a 3x3 kernel (correlation)."""
    p = np.pad(a, 1, mode="edge")
    out = np.zeros_like(a)
    for di in range(3):
        for dj in range(3):
            if k[di, dj] != 0:
                out = out + k[di, dj] * p[di:di + a.shape[0], dj:dj + a.shape[1]]
    return out


def edge_sobel(a):
    """MATLAB edge(a) with default arguments: method 'sobel', automatic threshold
    (cutoff = 4 * mean of the squared gradient magnitude), direction 'both', thinning on.

    MATLAB does the thinning in a compiled mex file. The rule below (non-maximum suppression
    along the dominant gradient direction, '<=' towards the left/upper neighbour and '<' towards
    the right/lower one) reproduces MATLAB's Q_FMI values exactly on bit-identical inputs
    (checked against values computed in MATLAB). Returns a logical map."""
    a = as_double(a)
    op = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]], dtype=np.float64) / 8
    bx = _imfilter_replicate(a, op.T)
    by = _imfilter_replicate(a, op)
    b = bx * bx + by * by
    cutoff = 4 * np.sum(b) / b.size
    bx, by = np.abs(bx), np.abs(by)
    m, n = b.shape
    left = b[:, np.r_[0, 0:n - 1]]
    right = b[:, np.r_[1:n, n - 1]]
    up = b[np.r_[0, 0:m - 1], :]
    down = b[np.r_[1:m, m - 1], :]
    horiz = (bx >= by) & (left <= b) & (right < b)
    vert = (by >= bx) & (up <= b) & (down < b)
    return (b > cutoff) & (horiz | vert)


# ---------------------------------------------------------------- window statistics
def _windows(x, w):
    """All (2w+1)x(2w+1) windows, flattened column-major like MATLAB aSub(:). Shape (N, L)."""
    k = 2 * w + 1
    v = sliding_window_view(x, (k, k))              # (m-2w, n-2w, k, k), row-major positions
    v = np.swapaxes(v, -1, -2)                      # column-major inside the window
    return v.reshape(-1, k * k)


def _normalise(win):
    mx, mn = win.max(axis=1, keepdims=True), win.min(axis=1, keepdims=True)
    const = (mx == mn)
    out = np.where(const, 1.0, (win - mn) / np.where(const, 1.0, mx - mn))
    return out


def _xlog2x(p):
    out = np.zeros_like(p)
    nz = p != 0
    out[nz] = -p[nz] * np.log2(np.abs(p[nz]))      # real(-p*log2(p)) also for p < 0
    return out


def _seqsum(x, axis=-1):
    """MATLAB-order (left-to-right) sum; numpy's pairwise sum changes the last bits, and the
    .m file branches on exact tests such as sum(aTemp.*fTemp) == 0."""
    return np.cumsum(x, axis=axis).take(-1, axis=axis)


def _sum2(win, k):
    """MATLAB sum(sum(X)) of a k x k window stored column-major: each column summed top to
    bottom, then the k column sums summed left to right."""
    cols = [_seqsum(win[:, j * k:(j + 1) * k]) for j in range(k)]
    tot = cols[0]
    for c in cols[1:]:
        tot = tot + c
    return tot


# order in which Q_FMI.m accumulates the joint entropy: (1,1); (i,1) i=2..l; (1,j) j=2..l;
# then i = 2..l (outer), j = 2..l (inner)
def _je_order(L):
    o = [(0, 0)] + [(i, 0) for i in range(1, L)] + [(0, j) for j in range(1, L)]
    o += [(i, j) for i in range(1, L) for j in range(1, L)]
    return np.array([i * L + j for i, j in o])


def _fmi_pairs(xn, fn):
    """fmi_xf for every window (rows), x and f already normalised; literal Q_FMI.m arithmetic."""
    L = xn.shape[1]
    k = int(round(np.sqrt(L)))
    xP = xn / _sum2(xn, k)[:, None]
    fP = fn / _sum2(fn, k)[:, None]
    xC = np.cumsum(xP, axis=1)
    fC = np.cumsum(fP, axis=1)
    xT = xP - (_seqsum(xP) / L)[:, None]
    fT = fP - (_seqsum(fP) / L)[:, None]
    s = _seqsum(xT * fT)
    den = np.sqrt(_seqsum(xT * xT) * _seqsum(fT * fT))
    c = np.where(s == 0, 0.0, s / np.where(s == 0, 1.0, den))
    i = np.arange(1, L + 1, dtype=np.float64)
    xV = _seqsum(i ** 2 * xP) - _seqsum(i * xP) ** 2
    fV = _seqsum(i ** 2 * fP) - _seqsum(i * fP) ** 2
    xSd = np.sqrt(np.maximum(xV, 0))
    fSd = np.sqrt(np.maximum(fV, 0))
    FI = fC[:, :, None]          # index i -> f
    XJ = xC[:, None, :]          # index j -> x
    up = c >= 0
    Mup = 0.5 * (FI + XJ - np.abs(FI - XJ))
    Mlo = 0.5 * (FI + XJ - 1 + np.abs(FI + XJ - 1))
    M = np.where(up[:, None, None], Mup, Mlo)
    cov = _seqsum((M - FI * XJ).reshape(M.shape[0], -1))       # i outer, j inner
    deg = (xSd == 0) | (fSd == 0) | (up & (c == 0))
    corr = cov / np.where(deg, 1.0, fSd * xSd)
    phi = np.where(deg, 0.0, c / np.where(corr == 0, 1.0, corr))
    Mp = np.zeros((M.shape[0], L + 1, L + 1))
    Mp[:, 1:, 1:] = M
    jUp = Mp[:, 1:, 1:] - Mp[:, :-1, 1:] - Mp[:, 1:, :-1] + Mp[:, :-1, :-1]
    # MATLAB writes the first-row/column terms with fewer subtractions; reproduce them exactly
    jUp[:, 0, 0] = M[:, 0, 0]
    jUp[:, 1:, 0] = M[:, 1:, 0] - M[:, :-1, 0]
    jUp[:, 0, 1:] = M[:, 0, 1:] - M[:, 0, :-1]
    jUp[:, 1:, 1:] = ((M[:, 1:, 1:] - M[:, :-1, 1:]) - M[:, 1:, :-1]) + M[:, :-1, :-1]
    jpdf = phi[:, None, None] * jUp + (1 - phi[:, None, None]) * (fP[:, :, None] * xP[:, None, :])
    JE = _seqsum(_xlog2x(jpdf).reshape(jpdf.shape[0], -1)[:, _je_order(L)])
    xE = _seqsum(_xlog2x(xP))
    fE = _seqsum(_xlog2x(fP))
    mi = (xE + fE) - JE
    return np.where(mi == 0, 0.0, 2 * mi / np.where(mi == 0, 1.0, xE + fE))


def Q_FMI(A, B, F, feature="edge", w=3, chunk=20000):
    """Feature mutual information (Q_FMI.m). Default feature 'edge' (MATLAB edge()), w = 3."""
    A, B, F = as_double(A), as_double(B), as_double(F)
    if feature == "edge":
        A, B, F = (edge_sobel(X).astype(np.float64) for X in (A, B, F))
    elif feature != "none":
        raise NotImplementedError("only feature='edge' (default) and 'none' are ported")
    w = w // 2
    aW, bW, fW = _windows(A, w), _windows(B, w), _windows(F, w)
    out = np.empty(aW.shape[0])
    for s in range(0, aW.shape[0], chunk):
        a, b, f = aW[s:s + chunk], bW[s:s + chunk], fW[s:s + chunk]
        eqA = np.all(a == f, axis=1)
        fmi_af = np.ones(a.shape[0])
        nA = ~eqA
        f_after = f.copy()
        if nA.any():
            fn = _normalise(f[nA])
            fmi_af[nA] = _fmi_pairs(_normalise(a[nA]), fn)
            f_after[nA] = fn                      # quirk: B-branch sees the normalised fSub
        eqB = np.all(b == f_after, axis=1)
        fmi_bf = np.ones(a.shape[0])
        nB = ~eqB
        if nB.any():
            fmi_bf[nB] = _fmi_pairs(_normalise(b[nB]), _normalise(f_after[nB]))
        out[s:s + chunk] = (fmi_af + fmi_bf) / 2
    return float(np.mean(out))


__all__ = ["Q_FMI", "edge_sobel"]
