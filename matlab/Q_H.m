function res = Q_H(im1, im2, fim)
% Q_H -- Hossny, Nahavandi & Creighton's Feature-Based Fusion Quality
% Metric ("FeatureQM", Algorithm 1 + eq. 13).
%
% Ref: M. Hossny, S. Nahavandi, D. Creighton, "Feature-Based Image
% Fusion Quality Metrics," ICIRA 2008, LNAI 5314, pp. 469-478 (see
% ../Python/978-3-540-88513-9_51_Q_H.pdf). The paper's own contribution
% over Piella's windowed UQI/SSIM (Sec. 2, eq. 10-12) is applying the
% quality metric Q0 to FEATURE maps s0(x), s0(f) of the source/fused
% images (e.g. local entropy) instead of to raw pixel intensities --
% "how much valuable [feature] information has not been transferred".
%
% Concrete choices (s0 = local entropy, Q0 = local SSIM, both 7x7
% windows) come from Z. Liu's metricHossny.m (July 2009), the only
% known reference implementation -- but that file does not actually
% run: it references undefined variables E_m1/E_m2 (should be
% E_im1/E_im2), and "E=ramda.*E1_(1-ramda).*E2" is a MATLAB syntax
% error (E1_ parsed as a call to an undefined function). Algorithm 1
% step 5 in the paper resolves both unambiguously:
%   e_{Q0,s}^w = ramda.*E1 + (1-ramda).*E2
% A third, more subtle bug: metricHossny.m casts im1/im2/fim to double
% BEFORE calling entropy()/nlfilter(...,'entropy') -- but MATLAB's
% entropy() treats a double-class array as if it already lies in
% [0,1] (imhist's default double behaviour), so on a normal [0,255]
% grayscale image every pixel falls in the single top histogram bin
% and the "entropy" map comes out ~0 everywhere (verified against real
% MATLAB: entropy(uint8_img) ~ 6.6, entropy(double(uint8_img)) ~ 0 on
% the same data). Fixed here by computing entropy on the image's
% original uint8 sample levels regardless of storage class.
%
% This function does NOT use nlfilter+wkeep (metricHossny.m's own
% approach: 'same'-sized filtering with zero-padded borders, then crop
% away the padding-affected border with wkeep) because wkeep's centering
% for a given odd/even crop amount is undocumented and error-prone to
% replicate in Python. Since windowSize=7 means only pixels within 3 of
% the border are ever padding-affected, and the reference code's own
% crop (bian=floor((7+1)/2)=4 pixels) already discards more than that,
% the padding convention never actually affects the final result --  so
% this implementation computes local entropy directly in "valid" mode
% (window fully inside the image, no padding at all), which is
% mathematically identical wherever it's defined and has an unambiguous,
% directly Python-portable size/centering convention.

windowSize = 7;   % s0 window (local entropy) AND Q0's SSIM Gaussian window
half = floor(windowSize / 2);  % = 3

im1 = double(im1);
im2 = double(im2);
fim = double(fim);

% ---- Step 1 (Algorithm 1): feature maps s0(x|w) = local entropy ----
% "valid" mode: output size (H-2*half, W-2*half), position (i,j) is the
% entropy of the windowSize x windowSize block centered there.
E_im1 = local_entropy_u8(im1, windowSize);
E_im2 = local_entropy_u8(im2, windowSize);
E_fim = local_entropy_u8(fim, windowSize);

% ---- Step 2: feature perceptual importance within each source image ----
p1 = E_im1 ./ max(E_im1(:));
p2 = E_im2 ./ max(E_im2(:));

% ---- Q0 = local SSIM applied to the ENTROPY feature maps (this
% paper's actual contribution vs. Piella's plain pixel-domain SSIM) ----
ssim_map1 = local_ssim_7x7(E_im1, E_fim);
ssim_map2 = local_ssim_7x7(E_im2, E_fim);
% filter2(...,'valid') with the 7x7 Gaussian shrinks by 2*half again ->
% crop p1/p2/E_im1/E_im2/E_fim (all currently the same size as
% E_im1/E_im2/E_fim above) down to ssim_map's smaller size by an equal
% margin on every side (unambiguous: the amount removed, 2*half, is
% always even for windowSize=7).
p1 = center_crop(p1, half);
p2 = center_crop(p2, half);
E_im1c = center_crop(E_im1, half);
E_im2c = center_crop(E_im2, half);
E_fimc = center_crop(E_fim, half);

% ---- Step 3: ex = px .* Q0(xs0, fs0), ey = py .* Q0(ys0, fs0) ----
E1 = p1 .* ssim_map1;
E2 = p2 .* ssim_map2;

% ---- Step 4: contribution/saliency weight lambda ----
% metricHossny.m's own concrete choice for the abstract saliency
% function s(x,f|w) = g(s0(x|w), s0(f|w)) of eq. 9 is the plain ratio
% g(a,b) = a/b (S_xf = E_im1./E_fim) -- computed here on the final,
% consistently-cropped entropy maps (metricHossny.m computes this
% BEFORE its own second crop, on mismatched sizes, which is why its
% "ramda.*E1" line could never actually run even after fixing the
% other two bugs).
eps_guard = 1e-12;  % avoid 0/0 in flat (zero-entropy) regions -- common
                     % in real medical images with large black background
S_xf = E_im1c ./ (E_fimc + eps_guard);
S_yf = E_im2c ./ (E_fimc + eps_guard);
ramda = S_xf ./ (S_xf + S_yf + eps_guard);

% ---- Step 5-6 + eq. 13: error map, then similarity = 1 - error ----
E = ramda .* E1 + (1 - ramda) .* E2;
res = mean2(1 - E);

end


function E = local_entropy_u8(img, windowSize)
% Local Shannon entropy, "valid" mode (no padding), 256-bin histogram
% over integer levels 0-255 -- matches real MATLAB's entropy() on a
% uint8 array (verified: nlfilter(uint8(img),[7 7],'entropy') equals
% this computed directly on the un-padded interior; entropy() on a
% double-class [0,255] array does NOT match this and must not be used).
lvl = uint8(min(max(round(img), 0), 255));
[H, W] = size(lvl);
half = floor(windowSize / 2);
outH = H - 2 * half;
outW = W - 2 * half;
E = zeros(outH, outW);
for i = 1:outH
    for j = 1:outW
        block = lvl(i:i + windowSize - 1, j:j + windowSize - 1);
        counts = accumarray(double(block(:)) + 1, 1, [256 1]);
        p = counts / sum(counts);
        p = p(p > 0);
        E(i, j) = -sum(p .* log2(p));
    end
end
end


function ssim_map = local_ssim_7x7(imgA, imgB)
% Same construction as metricHossny.m's local ssim_hossny() sub-function
% (NOT the shared ../../_deps_matlab/ssim_index.m, which defaults to an
% 11x11 window and the standard Wang-Bovik C1/C2): 7x7 Gaussian window,
% sigma=1.5, near-zero stabilization constants (this metric's own
% deliberate parameter choice, distinct from the standard SSIM
% defaults used elsewhere in this project).
window = fspecial('gaussian', 7, 1.5);
window = window / sum(window(:));
C1 = 2e-16;
C2 = 2e-16;

mu1 = filter2(window, imgA, 'valid');
mu2 = filter2(window, imgB, 'valid');
mu1_sq = mu1 .* mu1;
mu2_sq = mu2 .* mu2;
mu1_mu2 = mu1 .* mu2;
sigma1_sq = filter2(window, imgA .* imgA, 'valid') - mu1_sq;
sigma2_sq = filter2(window, imgB .* imgB, 'valid') - mu2_sq;
sigma12 = filter2(window, imgA .* imgB, 'valid') - mu1_mu2;

ssim_map = ((2 * mu1_mu2 + C1) .* (2 * sigma12 + C2)) ./ ...
    ((mu1_sq + mu2_sq + C1) .* (sigma1_sq + sigma2_sq + C2));
end


function out = center_crop(arr, margin)
[H, W] = size(arr);
out = arr(margin + 1:H - margin, margin + 1:W - margin);
end
