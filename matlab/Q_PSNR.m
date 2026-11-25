% Q_PSNR -- average Peak Signal-to-Noise Ratio: (PSNR(F,A)+PSNR(F,B))/2.
% NOT part of Liu et al. 2012's taxonomy, and doesn't fit any of its 4
% principle-based groups (pure pixelwise-MSE signal fidelity -- no
% information-theoretic, feature, structural, or HVS model involved),
% hence its own supplementary group here. Copied unchanged from
% Metrics_Matlab/Q_PSNR.m. Needs Image Processing Toolbox (psnr()).
function psnr_value = Q_PSNR(I_1, I_2, I_fused)

    psnr1 = psnr(I_fused,I_1,255);
    psnr2 = psnr(I_fused,I_2,255);

    psnr_value = (psnr1+psnr2)/2;

end
