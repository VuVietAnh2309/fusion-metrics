# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Port of Q_P.m (function Q_Zhao): Zhao et al.'s phase-congruency-based fusion metric,
with myphasecong3 (after P. Kovesi's phasecong3) and its lowpassfilter, as in the .m file.
"""
import numpy as np

from ._matlab import as_double, filter2_same, fspecial_gaussian


def _lowpassfilter_q_p(rows, cols, cutoff=0.45, n=15):
    """lowpassfilter inside Q_P.m (differs from the one in Q_FSIM.m: fix() grid and fftshift)."""
    x = (np.ones((rows, 1)) * np.arange(1, cols + 1) - (np.fix(cols / 2) + 1)) / cols
    y = (np.arange(1, rows + 1)[:, None] * np.ones((1, cols)) - (np.fix(rows / 2) + 1)) / rows
    radius = np.sqrt(x ** 2 + y ** 2)
    return np.fft.fftshift(1.0 / (1.0 + (radius / cutoff) ** (2 * n)))


def _range(n):
    if n % 2:
        return np.arange(-(n - 1) / 2, (n - 1) / 2 + 1) / (n - 1)
    return np.arange(-n / 2, n / 2) / n


def _myphasecong3(im, nscale=4, norient=6, minWaveLength=3, mult=2.1, sigmaOnf=0.55,
                  dThetaOnSigma=1.2, k=2.0, cutOff=0.5, g=10):
    """myphasecong3 in Q_P.m -> (phaseCongruency, M, m) (the orientation output is not used)."""
    im = as_double(im)
    epsilon = 0.0001
    thetaSigma = np.pi / norient / dThetaOnSigma
    rows, cols = im.shape
    imagefft = np.fft.fft2(im)
    zero = np.zeros((rows, cols))
    totalEnergy, totalSumAn = zero.copy(), zero.copy()
    covx2, covy2, covxy = zero.copy(), zero.copy(), zero.copy()
    x, y = np.meshgrid(_range(cols), _range(rows))
    radius = np.sqrt(x ** 2 + y ** 2)
    radius[rows // 2, cols // 2] = 1                  # radius(floor(rows/2)+1, floor(cols/2)+1) = 1
    theta = np.arctan2(-y, x)
    radius, theta = np.fft.ifftshift(radius), np.fft.ifftshift(theta)
    sintheta, costheta = np.sin(theta), np.cos(theta)
    lp = _lowpassfilter_q_p(rows, cols, 0.45, 15)
    logGabor = []
    for s in range(nscale):
        fo = 1.0 / (minWaveLength * mult ** s)
        lg = np.exp((-(np.log(radius / fo)) ** 2) / (2 * np.log(sigmaOnf) ** 2)) * lp
        lg[0, 0] = 0
        logGabor.append(lg)
    for o in range(norient):
        angl = o * np.pi / norient
        ds = sintheta * np.cos(angl) - costheta * np.sin(angl)
        dc = costheta * np.cos(angl) + sintheta * np.sin(angl)
        spread = np.exp((-np.abs(np.arctan2(ds, dc)) ** 2) / (2 * thetaSigma ** 2))
        sumE, sumO, sumAn, Energy = zero.copy(), zero.copy(), zero.copy(), zero.copy()
        EO, ifftF = [], []
        for s in range(nscale):
            filt = logGabor[s] * spread
            ifftF.append(np.real(np.fft.ifft2(filt)) * np.sqrt(rows * cols))
            eo = np.fft.ifft2(imagefft * filt)
            EO.append(eo)
            An = np.abs(eo)
            sumAn = sumAn + An
            sumE = sumE + np.real(eo)
            sumO = sumO + np.imag(eo)
            if s == 0:
                EM_n = np.sum(filt ** 2)
                maxAn = An
            else:
                maxAn = np.maximum(maxAn, An)
        XEnergy = np.sqrt(sumE ** 2 + sumO ** 2) + epsilon
        MeanE, MeanO = sumE / XEnergy, sumO / XEnergy
        for s in range(nscale):
            E, O = np.real(EO[s]), np.imag(EO[s])
            Energy = Energy + E * MeanE + O * MeanO - np.abs(E * MeanO - O * MeanE)
        medianE2n = np.median((np.abs(EO[0]) ** 2).T.ravel())
        meanE2n = -medianE2n / np.log(0.5)
        noisePower = meanE2n / EM_n
        EstSumAn2 = zero.copy()
        for s in range(nscale):
            EstSumAn2 = EstSumAn2 + ifftF[s] ** 2
        EstSumAiAj = zero.copy()
        for si in range(nscale - 1):
            for sj in range(si + 1, nscale):
                EstSumAiAj = EstSumAiAj + ifftF[si] * ifftF[sj]
        EstNoiseEnergy2 = 2 * noisePower * np.sum(EstSumAn2) + 4 * noisePower * np.sum(EstSumAiAj)
        tau = np.sqrt(EstNoiseEnergy2 / 2)
        EstNoiseEnergy = tau * np.sqrt(np.pi / 2)
        EstNoiseEnergySigma = np.sqrt((2 - np.pi / 2) * tau ** 2)
        T = (EstNoiseEnergy + k * EstNoiseEnergySigma) / 1.7
        Energy = np.maximum(Energy - T, zero)
        width = sumAn / (maxAn + epsilon) / nscale
        weight = 1.0 / (1 + np.exp((cutOff - width) * g))
        totalSumAn = totalSumAn + sumAn
        totalEnergy = totalEnergy + weight * Energy
        PC = weight * Energy / sumAn
        covx, covy = PC * np.cos(angl), PC * np.sin(angl)
        covx2 = covx2 + covx ** 2
        covy2 = covy2 + covy ** 2
        covxy = covxy + covx * covy
    phaseCongruency = totalEnergy / (totalSumAn + epsilon)
    covx2 = covx2 / (norient / 2)
    covy2 = covy2 / (norient / 2)
    covxy = covxy / norient
    denom = np.sqrt(covxy ** 2 + (covx2 - covy2) ** 2) + epsilon
    M = (covy2 + covx2 + denom) / 2
    m = (covy2 + covx2 - denom) / 2
    return phaseCongruency, M, m


def _correlation_coefficient(im1, im2, im_max, imf, mask1, mask2, mask3):
    """correlation_coeffcient in Q_P.m (11x11 Gaussian window, 'same' filtering, C = 1e-4)."""
    window = fspecial_gaussian(11, 1.5)
    window = window / window.sum()
    C = 0.0001
    im1, im2, im_max = mask1 * im1, mask2 * im2, mask3 * im_max
    f = lambda z: filter2_same(window, z)              # noqa: E731
    mu1, mu2, muf, mu_max = f(im1), f(im2), f(imf), f(im_max)
    sigma1_sq = f(im1 * im1) - mu1 * mu1
    sigma2_sq = f(im2 * im2) - mu2 * mu2
    sigmaMax_sq = f(im_max * im_max) - mu_max * mu_max
    sigmaf_sq = f(imf * imf) - muf * muf
    sigma1f = f(im1 * imf) - mu1 * muf
    sigma2f = f(im2 * imf) - mu2 * muf
    sigmaMaxf = f(im_max * imf) - mu_max * muf
    res = []
    for mask, sxf, sx in ((mask1, sigma1f, sigma1_sq), (mask2, sigma2f, sigma2_sq), (mask3, sigmaMaxf, sigmaMax_sq)):
        r = np.zeros_like(mu1)
        i = mask == 1
        r[i] = (sxf[i] + C) / (np.sqrt(np.abs(sx[i] * sigmaf_sq[i])) + C)
        res.append(r)
    result = np.maximum(np.maximum(res[0], res[1]), res[2])
    return float(np.sum(result) / np.sum(mask3))


def Q_P(A, B, F):
    """Zhao et al.'s phase-congruency metric (Q_P.m / Q_Zhao): product of the masked correlation
    scores of phase congruency and of the maximum and minimum moments, feature threshold 0.1."""
    fea_threshold = 0.1
    pc1, M1, m1 = _myphasecong3(A)
    pc2, M2, m2 = _myphasecong3(B)
    pcf, Mf, mf = _myphasecong3(F)
    mask = pc1 > pc2
    pc_max = mask * pc1 + (~mask) * pc2
    M_max = mask * M1 + (~mask) * M2
    m_max = mask * m1 + (~mask) * m2
    mask1 = (pc1 > fea_threshold).astype(np.float64)
    mask2 = (pc2 > fea_threshold).astype(np.float64)
    mask3 = (pc_max > fea_threshold).astype(np.float64)
    rPC = _correlation_coefficient(pc1, pc2, pc_max, pcf, mask1, mask2, mask3)
    rM = _correlation_coefficient(M1, M2, M_max, Mf, mask1, mask2, mask3)
    rm = _correlation_coefficient(m1, m2, m_max, mf, mask1, mask2, mask3)
    return float(rPC * rM * rm)


__all__ = ["Q_P"]
