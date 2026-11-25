# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Port of Q_FSIM.m: FSIM_TB(A, B, F) = (FSIM(F, A) + FSIM(F, B)) / 2, grayscale path,
with phasecong2 and lowpassfilter exactly as in the .m file."""
import numpy as np

from ._matlab import as_double, conv2_same


def _grid(rows, cols):
    if cols % 2:
        xr = np.arange(-(cols - 1) / 2, (cols - 1) / 2 + 1) / (cols - 1)
    else:
        xr = np.arange(-cols / 2, cols / 2) / cols
    if rows % 2:
        yr = np.arange(-(rows - 1) / 2, (rows - 1) / 2 + 1) / (rows - 1)
    else:
        yr = np.arange(-rows / 2, rows / 2) / rows
    return np.meshgrid(xr, yr)


def _lowpassfilter(rows, cols, cutoff, n):
    x, y = _grid(rows, cols)
    radius = np.sqrt(x ** 2 + y ** 2)
    return np.fft.ifftshift(1.0 / (1.0 + (radius / cutoff) ** (2 * n)))


def _phasecong2(im):
    nscale, norient = 4, 4
    minWaveLength, mult, sigmaOnf = 6, 2, 0.55
    dThetaOnSigma, k, epsilon = 1.2, 2.0, 0.0001
    thetaSigma = np.pi / norient / dThetaOnSigma
    rows, cols = im.shape
    imagefft = np.fft.fft2(im)
    x, y = _grid(rows, cols)
    radius = np.fft.ifftshift(np.sqrt(x ** 2 + y ** 2))
    theta = np.fft.ifftshift(np.arctan2(-y, x))
    radius[0, 0] = 1
    sintheta, costheta = np.sin(theta), np.cos(theta)
    lp = _lowpassfilter(rows, cols, 0.45, 15)
    logGabor = []
    for s in range(nscale):
        fo = 1.0 / (minWaveLength * mult ** s)
        lg = np.exp((-(np.log(radius / fo)) ** 2) / (2 * np.log(sigmaOnf) ** 2)) * lp
        lg[0, 0] = 0
        logGabor.append(lg)
    spread = []
    for o in range(norient):
        angl = o * np.pi / norient
        ds = sintheta * np.cos(angl) - costheta * np.sin(angl)
        dc = costheta * np.cos(angl) + sintheta * np.sin(angl)
        dtheta = np.abs(np.arctan2(ds, dc))
        spread.append(np.exp((-dtheta ** 2) / (2 * thetaSigma ** 2)))
    EnergyAll = np.zeros((rows, cols))
    AnAll = np.zeros((rows, cols))
    for o in range(norient):
        sumE = np.zeros((rows, cols)); sumO = np.zeros((rows, cols)); sumAn = np.zeros((rows, cols))
        Energy = np.zeros((rows, cols))
        EO, ifftF = [], []
        for s in range(nscale):
            filt = logGabor[s] * spread[o]
            ifftF.append(np.real(np.fft.ifft2(filt)) * np.sqrt(rows * cols))
            eo = np.fft.ifft2(imagefft * filt)
            EO.append(eo)
            An = np.abs(eo)
            sumAn = sumAn + An
            sumE = sumE + np.real(eo)
            sumO = sumO + np.imag(eo)
            if s == 0:
                EM_n = np.sum(filt ** 2)
        XEnergy = np.sqrt(sumE ** 2 + sumO ** 2) + epsilon
        MeanE = sumE / XEnergy
        MeanO = sumO / XEnergy
        for s in range(nscale):
            E, O = np.real(EO[s]), np.imag(EO[s])
            Energy = Energy + E * MeanE + O * MeanO - np.abs(E * MeanO - O * MeanE)
        medianE2n = np.median((np.abs(EO[0]) ** 2).T.ravel())
        meanE2n = -medianE2n / np.log(0.5)
        noisePower = meanE2n / EM_n
        EstSumAn2 = np.zeros((rows, cols))
        for s in range(nscale):
            EstSumAn2 = EstSumAn2 + ifftF[s] ** 2
        EstSumAiAj = np.zeros((rows, cols))
        for si in range(nscale - 1):
            for sj in range(si + 1, nscale):
                EstSumAiAj = EstSumAiAj + ifftF[si] * ifftF[sj]
        EstNoiseEnergy2 = 2 * noisePower * np.sum(EstSumAn2) + 4 * noisePower * np.sum(EstSumAiAj)
        tau = np.sqrt(EstNoiseEnergy2 / 2)
        EstNoiseEnergy = tau * np.sqrt(np.pi / 2)
        EstNoiseEnergySigma = np.sqrt((2 - np.pi / 2) * tau ** 2)
        T = (EstNoiseEnergy + k * EstNoiseEnergySigma) / 1.7
        Energy = np.maximum(Energy - T, 0)
        EnergyAll = EnergyAll + Energy
        AnAll = AnAll + sumAn
    return EnergyAll / AnAll


def _fsim(ref, dis):
    Y1, Y2 = as_double(ref), as_double(dis)
    rows, cols = Y1.shape
    Fd = max(1, int(np.floor(min(rows, cols) / 256 + 0.5)))
    if Fd > 1:
        k = np.ones((Fd, Fd)) / (Fd * Fd)
        Y1 = conv2_same(Y1, k)[::Fd, ::Fd]
        Y2 = conv2_same(Y2, k)[::Fd, ::Fd]
    PC1, PC2 = _phasecong2(Y1), _phasecong2(Y2)
    dx = np.array([[3, 0, -3], [10, 0, -10], [3, 0, -3]]) / 16
    dy = np.array([[3, 10, 3], [0, 0, 0], [-3, -10, -3]]) / 16
    g1 = np.sqrt(conv2_same(Y1, dx) ** 2 + conv2_same(Y1, dy) ** 2)
    g2 = np.sqrt(conv2_same(Y2, dx) ** 2 + conv2_same(Y2, dy) ** 2)
    T1, T2 = 0.85, 160
    PCSim = (2 * PC1 * PC2 + T1) / (PC1 ** 2 + PC2 ** 2 + T1)
    GSim = (2 * g1 * g2 + T2) / (g1 ** 2 + g2 ** 2 + T2)
    PCm = np.maximum(PC1, PC2)
    return np.sum(GSim * PCSim * PCm) / np.sum(PCm)


def Q_FSIM(A, B, F):
    """Feature similarity for fusion (Q_FSIM.m / FSIM_TB): mean of FSIM(F, A) and FSIM(F, B)."""
    return float((_fsim(F, A) + _fsim(F, B)) / 2)


__all__ = ["Q_FSIM"]
