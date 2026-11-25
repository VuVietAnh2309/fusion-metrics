function res=Q_C(im1,im2,fim)
% Q_C -- Liu et al. 2012's "Cvejie's Metric" (Sec. 2.3.2, eq. 44-45).
%
% Ref: N. Cvejic, A. Loza, D. Bull, and N. Canagarajah, "A Similarity
% Metric for Assessment of Image Fusion Algorithms," Int'l J. Signal
% Processing, vol. 2, no. 3, pp. 178-182, 2005.
%
% Faithful port of metricCvejic.m's "case 2" branch (Z. Liu @NRCC, July
% 2009) -- the reference implementation's "case 1" is an unfinished
% stub in the original file (hardcoded res=0, marked "to be finished"),
% so only case 2 (which matches this paper's eq. 44-45 exactly) is
% used. Renamed to match this project's Q_<NAME>(A,B,F) convention.
%
% Uses a local UIQI-like windowed similarity (Gaussian window, 7x7,
% near-zero stabilizing constants C1=C2=2e-16) for Q(A,F|w), Q(B,F|w).
%
% Zero-denominator guard (not in the original metricCvejic.m): wherever
% A and B are both locally flat (e.g. black background in real medical
% images), sigma_XF+sigma_YF is exactly 0, making simXYF = 0/0 = NaN,
% which poisons the whole-image mean2(). Fixed the same way this
% project's sibling metric Q_Peilla.m already handles the identical
% pattern: add 0.5 to both sigmas when their sum is 0, so simXYF=0.5
% (even blend -- neither source has local texture to prefer there).

im1=double(im1);
im2=double(im2);
fim=double(fim);

[mssim2, ssim_map2, sigma_XF] = ssim_yang(im1, fim); %#ok<ASGLU>
[mssim3, ssim_map3, sigma_YF] = ssim_yang(im2, fim); %#ok<ASGLU>

buffer=sigma_XF+sigma_YF;
test=(buffer==0); test=test*0.5;
sigma_XF=sigma_XF+test; sigma_YF=sigma_YF+test;

simXYF=sigma_XF./(sigma_XF+sigma_YF);
sim=simXYF.*ssim_map2+(1-simXYF).*ssim_map3;

index=find(simXYF<0);
sim(index)=0;

index=find(simXYF>1);
sim(index)=1;
res=mean2(sim);

end

%%
function [mssim, ssim_map, sigma12] = ssim_yang(img1, img2)

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
