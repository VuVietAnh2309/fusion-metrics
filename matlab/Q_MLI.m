% Q_MLI -- Mean Local Intensity. NOT part of Liu et al. 2012's 12-metric
% taxonomy (single-image, no-reference statistic).
% Copied unchanged from Metrics_Matlab/Q_MLI.m.
function MLI = Q_MLI(F)
% Q_CI - Contrast Index

F = double(F);

if min(F(:)) >= 0 && max(F(:)) <= 1
    % Image in [0,1]
    MLI = mean(F(:));

elseif min(F(:)) >= 0 && max(F(:)) <= 255
    % Image in [0,255]
    MLI = mean(F(:) / 255);

else
    error('Input image must have values in either [0,1] or [0,255].');
end

end