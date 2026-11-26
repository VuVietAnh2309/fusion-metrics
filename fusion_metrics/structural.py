# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Ports of MS_SSIM.m (SSIM, Z. Wang 2003) and Q_UIQI.m (UIQI, Z. Wang & A. C. Bovik 2001).
Argument order (A, B, F). Both are the mean of the index between each source and F.
"""
import numpy as np

from ._matlab import as_double, filter2_valid, fspecial_gaussian, mtimes_elem


# ---------------------------------------------------------------- MS_SSIM.m (ssim_index)
def _ssim_index(img1, img2, K=(0.01, 0.03), L=255.0):
    """ssim_index.m with its default arguments: 11x11 Gaussian window (sigma 1.5), K=[0.01 0.03],
    L=255, 'valid' filtering; returns mean2(ssim_map)."""
    M, N = np.shape(img1)
    if M < 11 or N < 11:
        # the .m file returns without assigning its output, which is an error in MATLAB
        raise ValueError("MS_SSIM: images must be at least 11x11")
    window = fspecial_gaussian(11, 1.5)
    window = window / window.sum()
    C1, C2 = (K[0] * L) ** 2, (K[1] * L) ** 2
    img1, img2 = as_double(img1), as_double(img2)
    mu1 = filter2_valid(window, img1)
    mu2 = filter2_valid(window, img2)
    mu1_sq, mu2_sq, mu1_mu2 = mu1 * mu1, mu2 * mu2, mu1 * mu2
    sigma1_sq = filter2_valid(window, img1 * img1) - mu1_sq
    sigma2_sq = filter2_valid(window, img2 * img2) - mu2_sq
    sigma12 = filter2_valid(window, img1 * img2) - mu1_mu2
    ssim_map = ((2 * mu1_mu2 + C1) * (2 * sigma12 + C2)) / ((mu1_sq + mu2_sq + C1) * (sigma1_sq + sigma2_sq + C2))
    return float(np.mean(ssim_map))


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


__all__ = ["MS_SSIM", "Q_UIQI"]
