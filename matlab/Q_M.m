function res = Q_M(im1, im2, fim)
% Q_M -- Liu et al. 2012's Multiscale Metric (Sec. 2.2.2, eq. 22-25).
%
% Ref: Wang and Liu's two-level Haar wavelet fusion metric, as
% described in Liu et al. 2012 [25]. Written from scratch from the
% paper's equations (no reference toolbox implementation was found --
% the codebase's "metricWang.m" is a mislabeled NCIE implementation,
% not this metric).
%
% Uses MATLAB's own dwt2(...,'haar') (Wavelet Toolbox, default 'sym'
% symmetric-extension boundary mode) for the two-level decomposition,
% so this *is* the ground-truth reference for the Python port, not a
% copy of a third-party author's code.
%
% Per-level weights alpha_s (eq. 16 of the original source paper --
% Peng-wei Wang & Bo Liu, "A Novel Image Fusion Metric Based on
% Multi-Scale Analysis," ICSP2008 -- = eq. 25 of Liu et al. 2012):
% Section IV of that paper states these explicitly -- "We empirically
% set the number of wavelet decomposition levels to 2, i.e., N = 2,
% and the constants alpha_1 and alpha_2 ... are empirically set to 2/3
% and 1/3, respectively." (Liu et al. 2012's own summarizing text is
% vaguer -- "alpha_s is a constant to adjust the relative importance
% of different scales" -- which led an earlier version of this file to
% wrongly assume the value was unspecified and default to equal
% weighting; corrected 2026-08-28 after checking the original source
% paper directly.)

N = 2;
alpha_s = [2/3, 1/3];

im1 = double(im1);
im2 = double(im2);
fim = double(fim);

Q_levels = zeros(1, N);
LL1 = im1; LL2 = im2; LLf = fim;
for s = 1:N
    [LL1, LH1, HL1, HH1] = dwt2(LL1, 'haar');
    [LL2, LH2, HL2, HH2] = dwt2(LL2, 'haar');
    [LLf, LHf, HLf, HHf] = dwt2(LLf, 'haar');

    psiH_AF = exp(-abs(LH1 - LHf));
    psiV_AF = exp(-abs(HL1 - HLf));
    psiD_AF = exp(-abs(HH1 - HHf));
    EP_AF = (psiH_AF + psiV_AF + psiD_AF) / 3;

    psiH_BF = exp(-abs(LH2 - LHf));
    psiV_BF = exp(-abs(HL2 - HLf));
    psiD_BF = exp(-abs(HH2 - HHf));
    EP_BF = (psiH_BF + psiV_BF + psiD_BF) / 3;

    wA = LH1.^2 + HL1.^2 + HH1.^2;
    wB = LH2.^2 + HL2.^2 + HH2.^2;

    num = sum(sum(EP_AF .* wA + EP_BF .* wB));
    den = sum(sum(wA + wB));
    Q_levels(s) = num / den;
end

res = prod(Q_levels .^ alpha_s);

end
