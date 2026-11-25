% Q_CC -- average Pearson correlation coefficient between the fused
% image and each source: (corr(A,F)+corr(B,F))/2. NOT part of Liu et al.
% 2012's taxonomy, and doesn't fit any of its 4 principle-based groups
% (pure pixel-level linear-correlation fidelity -- no
% information-theoretic, feature, structural, or HVS model involved),
% same rationale as this group's Q_PSNR. User-provided; copied from
% Metrics_Matlab/Q_CC.m with one bug fix: the original declared its
% output variable as Q_CC (function Q_CC = metric_CC(...)) but the body
% computed into a differently-named local variable CC and never assigned
% Q_CC, so the function threw "Output argument 'Q_CC' ... not assigned a
% value" and never actually ran. Fixed by assigning the result straight
% into the declared output variable.
function Q_CC = metric_CC(A, B, F)

    A = double(A(:));
    B = double(B(:));
    F = double(F(:));

    C1 = corrcoef(A, F);
    C2 = corrcoef(B, F);

    r_AF = C1(1,2);
    r_BF = C2(1,2);

    Q_CC = (r_AF + r_BF) / 2;

end
