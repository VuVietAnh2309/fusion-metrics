% Q_AG -- Average Gradient (sharpness). NOT part of Liu et al. 2012's
% 12-metric taxonomy (single-image, no-reference statistic -- takes only
% the fused image F, not the (A,B,F) triple the paper's metrics compare).
% Copied unchanged from Metrics_Matlab/Q_AG.m.
function AG = Q_AG(F)
% Q_CI - Contrast Index

F = double(F);

if min(F(:)) >= 0 && max(F(:)) <= 1
    % Image in [0,1]
    AG = sharpness(F);

elseif min(F(:)) >= 0 && max(F(:)) <= 255
    % Image in [0,255]
    AG = sharpness(F./255);

else
    error('Input image must have values in either [0,1] or [0,255].');
end

end

function f = sharpness(I)

[Gx, Gy] = imgradientxy(I,'intermediate');

S=sqrt(Gx.*Gx+Gy.*Gy);
f = sum(sum(S))./(numel(Gx));

end