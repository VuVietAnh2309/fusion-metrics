function [CE, D1, D2] = Q_CE(A, B, F)
%CALCULATE_CE Calculate histogram-based Cross Entropy metric
%
%   CE = 0.5 * [ D_KL(hA || hF) + D_KL(hB || hF) ]
%
% where:
%
%   D_KL(P || Q) = sum P(i) * log2(P(i)/Q(i))
%
% Inputs:
%   A : first source image
%   B : second source image
%   F : fused image
%
% Outputs:
%   CE : final Cross Entropy metric
%   D1 : D_KL(A || F)
%   D2 : D_KL(B || F)
%
% Lower CE is better.


    %% Convert images to grayscale uint8
    A = prepare_image(A);
    B = prepare_image(B);
    F = prepare_image(F);


    %% Compute normalized histograms
    hA = imhist(A, 256);
    hB = imhist(B, 256);
    hF = imhist(F, 256);

    hA = double(hA);
    hB = double(hB);
    hF = double(hF);

    hA = hA / sum(hA);
    hB = hB / sum(hB);
    hF = hF / sum(hF);


    %% Numerical stability
    epsVal = 1e-12;


    %% KL divergence A -> F
    maskA = hA > 0;

    D1 = sum( ...
        hA(maskA) .* ...
        log2( ...
            hA(maskA) ./ ...
            max(hF(maskA), epsVal) ...
        ) ...
    );


    %% KL divergence B -> F
    maskB = hB > 0;

    D2 = sum( ...
        hB(maskB) .* ...
        log2( ...
            hB(maskB) ./ ...
            max(hF(maskB), epsVal) ...
        ) ...
    );


    %% Final CE
    CE = 0.5 * (D1 + D2);

end


function img = prepare_image(img)
%PREPARE_IMAGE Convert image to grayscale uint8

    %% RGB to grayscale
    if ndims(img) == 3
        img = rgb2gray(img);
    end

    %% Convert depending on data type/range
    if isa(img, 'uint8')
        return;
    end

    img = double(img);

    %% Normalized image [0,1]
    if max(img(:)) <= 1
        img = img * 255;
    end

    %% Clip
    img(img < 0) = 0;
    img(img > 255) = 255;

    img = uint8(img);

end