# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Small helpers that reproduce MATLAB built-ins exactly (as used by the .m metrics)."""
import cv2
import numpy as np
from scipy.signal import convolve2d, correlate2d

EPS = np.finfo(np.float64).eps          # MATLAB eps
RGB2GRAY = (0.298936021293775, 0.587043074451121, 0.114020904255103)  # MATLAB rgb2gray


def mround(x):
    """MATLAB round: halves away from zero (numpy rounds halves to even)."""
    x = np.asarray(x, dtype=np.float64)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def to_uint8(x):
    """MATLAB uint8(x): round half away from zero, saturate to [0, 255]."""
    return np.clip(mround(x), 0, 255).astype(np.uint8)


def rgb2gray(rgb):
    """MATLAB rgb2gray for a uint8 RGB array (H, W, 3) in R, G, B order -> uint8."""
    rgb = rgb.astype(np.float64)
    y = RGB2GRAY[0] * rgb[..., 0] + RGB2GRAY[1] * rgb[..., 1] + RGB2GRAY[2] * rgb[..., 2]
    return to_uint8(y)


def read_gray(path):
    """MATLAB imread + rgb2gray, the input convention of the reference .m implementations.

    A single-channel file is returned unchanged; a 3-channel file (pseudo-colour PET/SPECT or a
    colour fused image) is converted with MATLAB rgb2gray. Returns uint8 (H, W).
    """
    im = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if im is None:
        raise FileNotFoundError(path)
    if im.ndim == 2:
        return im.astype(np.uint8)
    if im.shape[2] == 4:
        im = im[..., :3]
    return rgb2gray(im[..., ::-1])  # OpenCV is BGR


def read_gray_cv(path):
    """OpenCV convention common in Python fusion code: cv2 BGR2GRAY (0.299, 0.587, 0.114,
    rounded) scaled to [0,1], then x*255 truncated to uint8. Use it to match numbers produced by
    such code; use read_gray() to match the MATLAB reference implementations."""
    im = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if im is None:
        raise FileNotFoundError(path)
    if im.ndim == 3:
        im = cv2.cvtColor(im[..., :3], cv2.COLOR_BGR2GRAY)
    x = im.astype(np.float64) / 255.0
    return np.clip(x * 255.0, 0, 255).astype(np.uint8)


def conv2_same(x, k):
    """MATLAB conv2(x, k, 'same'): the central part of the full convolution, starting at index
    floor(size(k)/2) in each dimension (scipy's mode='same' differs by one pixel for even kernels)."""
    k = np.asarray(k, dtype=np.float64)
    full = convolve2d(x, k, mode="full")
    r0, c0 = k.shape[0] // 2, k.shape[1] // 2
    return full[r0:r0 + x.shape[0], c0:c0 + x.shape[1]]


def filter2_valid(k, x):
    """MATLAB filter2(k, x, 'valid') (correlation)."""
    return correlate2d(x, k, mode="valid")


def fspecial_gaussian(n, sigma):
    """MATLAB fspecial('gaussian', n, sigma) for a square odd/even size n."""
    siz = (n - 1) / 2.0
    y, x = np.meshgrid(np.arange(-siz, siz + 1), np.arange(-siz, siz + 1), indexing="ij")
    h = np.exp(-(x * x + y * y) / (2.0 * sigma * sigma))
    h[h < EPS * h.max()] = 0
    s = h.sum()
    return h / s if s != 0 else h


def as_double(x):
    return np.asarray(x).astype(np.float64)
