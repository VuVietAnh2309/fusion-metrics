
function NABF = Q_NABF(A, B, F)
%CALCULATE_NABF
% Modified fusion-artifact measure N_AB/F
%
% Based on:
% B. K. Shreyamsha Kumar,
% "Multifocus and Multispectral Image Fusion based on Pixel
% Significance using Discrete Cosine Harmonic Wavelet Transform"
%
% Lower NABF indicates fewer fusion-generated artifacts.
%
% Inputs:
%   A : source image 1
%   B : source image 2
%   F : fused image
%
% Output:
%   NABF : modified fusion artifact measure


%% ------------------------------------------------------------
% Parameters used in Kumar/Petrovic implementation
%% ------------------------------------------------------------

Td = 2;

wt_min = 0.001;

Lg = 1.5;

Nrg = 0.9999;
kg = 19;
sigmag = 0.5;

Nra = 0.9995;
ka = 22;
sigmaa = 0.5;


%% ------------------------------------------------------------
% Prepare images
%% ------------------------------------------------------------

A = prepare_gray(A);
B = prepare_gray(B);
F = prepare_gray(F);

if ~isequal(size(A), size(B), size(F))
    error('A, B and F must have identical dimensions.');
end


%% ------------------------------------------------------------
% Sobel edge components
%% ------------------------------------------------------------

[gvA, ghA] = sobel_kumar(A);
[gvB, ghB] = sobel_kumar(B);
[gvF, ghF] = sobel_kumar(F);


%% Edge strength

gA = sqrt(ghA.^2 + gvA.^2);
gB = sqrt(ghB.^2 + gvB.^2);
gF = sqrt(ghF.^2 + gvF.^2);


%% ------------------------------------------------------------
% Relative edge strengths
%% ------------------------------------------------------------

gAF = zeros(size(A));
gBF = zeros(size(B));

maskAF = (gA > 0) & (gF > 0);

gAF(maskAF) = ...
    min(gA(maskAF), gF(maskAF)) ./ ...
    max(gA(maskAF), gF(maskAF));

maskBF = (gB > 0) & (gF > 0);

gBF(maskBF) = ...
    min(gB(maskBF), gF(maskBF)) ./ ...
    max(gB(maskBF), gF(maskBF));


%% ------------------------------------------------------------
% Edge orientations
%
% IMPORTANT:
% Original implementation uses atan(gv/gh),
% not atan2().
%% ------------------------------------------------------------

aA = zeros(size(A));
aB = zeros(size(B));
aF = zeros(size(F));

maskA = ~((gvA == 0) & (ghA == 0));
maskB = ~((gvB == 0) & (ghB == 0));
maskF = ~((gvF == 0) & (ghF == 0));

aA(maskA) = atan( ...
    gvA(maskA) ./ ghA(maskA) ...
);

aB(maskB) = atan( ...
    gvB(maskB) ./ ghB(maskB) ...
);

aF(maskF) = atan( ...
    gvF(maskF) ./ ghF(maskF) ...
);


%% ------------------------------------------------------------
% Relative edge orientation
%% ------------------------------------------------------------

aAF = abs( ...
    abs(aA - aF) - pi/2 ...
) * 2/pi;

aBF = abs( ...
    abs(aB - aF) - pi/2 ...
) * 2/pi;


%% ------------------------------------------------------------
% Edge-strength preservation
%% ------------------------------------------------------------

QgAF = Nrg ./ ...
    (1 + exp(-kg * (gAF - sigmag)));

QgBF = Nrg ./ ...
    (1 + exp(-kg * (gBF - sigmag)));


%% ------------------------------------------------------------
% Edge-orientation preservation
%% ------------------------------------------------------------

QaAF = Nra ./ ...
    (1 + exp(-ka * (aAF - sigmaa)));

QaBF = Nra ./ ...
    (1 + exp(-ka * (aBF - sigmaa)));


%% ------------------------------------------------------------
% Total edge preservation
%% ------------------------------------------------------------

QAF = sqrt(QgAF .* QaAF);

QBF = sqrt(QgBF .* QaBF);


%% ------------------------------------------------------------
% Perceptual weights
%% ------------------------------------------------------------

wtA = wt_min * ones(size(A));
wtB = wt_min * ones(size(B));

maskWA = gA >= Td;

wtA(maskWA) = ...
    gA(maskWA).^Lg;

maskWB = gB >= Td;

wtB(maskWB) = ...
    gB(maskWB).^Lg;


%% Total source-edge weight

wt_sum = sum(wtA(:) + wtB(:));

if wt_sum == 0
    NABF = 0;
    return;
end


%% ------------------------------------------------------------
% Artifact mask
%
% AM = 1 where fused edge is stronger than BOTH source edges.
%% ------------------------------------------------------------

AM = double( ...
    (gF > gA) & ...
    (gF > gB) ...
);


%% ------------------------------------------------------------
% Modified N_AB/F
%% ------------------------------------------------------------

artifactMap = ...
    AM .* ...
    ( ...
        (1 - QAF) .* wtA + ...
        (1 - QBF) .* wtB ...
    );

NABF = sum(artifactMap(:)) / wt_sum;

end


% ==============================================================
% Image preparation
% ==============================================================

function img = prepare_gray(img)

    if ndims(img) == 3
        img = rgb2gray(img);
    end

    img = double(img);

end


% ==============================================================
% Sobel implementation used by Kumar
% ==============================================================

function [gv, gh] = sobel_kumar(x)

    vtemp = [ ...
        -1  0  1;
        -2  0  2;
        -1  0  1 ...
    ] / 8;

    htemp = [ ...
        -1 -2 -1;
         0  0  0;
         1  2  1 ...
    ] / 8;

    gv = imfilter( ...
        x, ...
        vtemp, ...
        'symmetric', ...
        'same', ...
        'corr' ...
    );

    gh = imfilter( ...
        x, ...
        htemp, ...
        'symmetric', ...
        'same', ...
        'corr' ...
    );

end