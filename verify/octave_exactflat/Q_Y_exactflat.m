% VERIFICATION ONLY -- copy of ../../matlab/Q_Y.m with one change: in ssim_yang, windows in
% which an image is constant get exactly zero variance/covariance. With C1 = C2 = 2e-16 the
% original file returns round-off-dependent values there (MATLAB, Octave and NumPy disagree);
% this copy is the reference for the Python port, which applies the same rule.
function res=Q_Y_exactflat(im1,im2,fim)
% Q_Y -- Liu et al. 2012's "Yang's Metric" (Sec. 2.3.3, eq. 46).
%
% Ref: C. Yang, J. Zhang, X. Wang, and X. Liu, "A novel similarity
% based quality metric for image fusion," Information Fusion, vol. 9,
% pp. 156-160, 2008.
%
% Faithful port of metricYang.m (Z. Liu @NRCC, July 2009). Renamed to
% match this project's Q_<NAME>(A,B,F) convention.
%
% Zero-denominator guard (not in the original metricYang.m): wherever
% A and B are both locally flat (e.g. black background in real medical
% images), sigma1_sq1+sigma2_sq1 is exactly 0, making ramda = 0/0 = NaN,
% which poisons the whole-image mean2(). Fixed the same way this
% project's sibling metric Q_Peilla.m already handles the identical
% pattern: add 0.5 to both sigmas when their sum is 0, so ramda=0.5
% (even blend -- neither source has local texture to prefer there).

im1=double(im1);
im2=double(im2);
fim=double(fim);

[mssim1, ssim_map1, sigma1_sq1,sigma2_sq1] = ssim_yang(im1, im2); %#ok<ASGLU>
[mssim2, ssim_map2, sigma1_sq2,sigma2_sq2] = ssim_yang(im1, fim); %#ok<ASGLU>
[mssim3, ssim_map3, sigma1_sq3,sigma2_sq3] = ssim_yang(im2, fim); %#ok<ASGLU>

bin_map=ssim_map1>=0.75;

buffer=sigma1_sq1+sigma2_sq1;
test=(buffer==0); test=test*0.5;
sigma1_sq1=sigma1_sq1+test; sigma2_sq1=sigma2_sq1+test;

ramda=sigma1_sq1./(sigma1_sq1+sigma2_sq1);

Q1=(ramda.*ssim_map2+(1-ramda).*ssim_map3).*bin_map;

Q2=(max(ssim_map2,ssim_map3)).*(~bin_map);

Q=mean2(Q1+Q2);

res=Q;

end

%%
function [mssim, ssim_map, sigma1_sq,sigma2_sq] = ssim_yang(img1, img2)

[M N] = size(img1); %#ok<NASGU>
if ((M < 11) | (N < 11))
   ssim_map = -Inf;
  return
end
window = fspecial('gaussian', 7, 1.5);

C1 = 2e-16;
C2 = 2e-16;

window = window/sum(sum(window));
img1 = double(img1);
img2 = double(img2);
mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;
% --- verification-only change (not in the reference file): windows in which an image is
% constant are evaluated in exact arithmetic (variance and covariance = 0), as in the Python port.
f1 = flatwin(img1, 7); f2 = flatwin(img2, 7);
sigma1_sq(f1) = 0; sigma2_sq(f2) = 0; sigma12(f1 | f2) = 0;
if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
    denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));

   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

end

function f = flatwin(img, w)
C = im2col(img, [w w], 'sliding');
f = reshape(max(C, [], 1) == min(C, [], 1), size(img, 1) - w + 1, size(img, 2) - w + 1);
end
