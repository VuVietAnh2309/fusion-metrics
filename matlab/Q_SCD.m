
function [SCD, SCD_A, SCD_B] = Q_SCD(A, B, F)
%CALCULATE_SCD Calculate Sum of Correlations of Differences
%
% Formula:
%
%   SCD = CC(F-B, A) + CC(F-A, B)
%
% Inputs:
%   A : source image 1
%   B : source image 2
%   F : fused image
%
% Outputs:
%   SCD   : final SCD value
%   SCD_A : CC(F-B, A)
%   SCD_B : CC(F-A, B)
%
% Higher SCD is better.


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
% Difference images
%
% IMPORTANT:
% Images are already converted to double,
% so subtraction does not cause uint8 underflow.
%% ------------------------------------------------------------

diffA = F - B;
diffB = F - A;


%% ------------------------------------------------------------
% Correlations
%% ------------------------------------------------------------

SCD_A = correlation_coefficient(diffA, A);

SCD_B = correlation_coefficient(diffB, B);


%% ------------------------------------------------------------
% Final SCD
%% ------------------------------------------------------------

SCD = SCD_A + SCD_B;

end


% ==============================================================
% Pearson correlation coefficient
% ==============================================================

function CC = correlation_coefficient(X, Y)

    X = double(X);
    Y = double(Y);

    if ~isequal(size(X), size(Y))
        error('X and Y must have identical dimensions.');
    end

    meanX = mean(X(:));
    meanY = mean(Y(:));

    Xc = X - meanX;
    Yc = Y - meanY;

    numerator = sum(Xc(:) .* Yc(:));

    denominator = sqrt( ...
        sum(Xc(:).^2) * ...
        sum(Yc(:).^2) ...
    );

    if denominator < eps
        CC = 0;
        return;
    end

    CC = numerator / denominator;

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