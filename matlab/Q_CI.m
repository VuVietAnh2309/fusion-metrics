% Q_CI -- Contrast Index (intensity std. dev.). NOT part of Liu et al.
% 2012's 12-metric taxonomy (single-image, no-reference statistic).
% Copied unchanged from Metrics_Matlab/Q_CI.m.
function Q_CI = Q_CI(F)
% Q_CI - Contrast Index

F = double(F);

if min(F(:)) >= 0 && max(F(:)) <= 1
    % Image in [0,1]
    Q_CI = std(F(:), 1);

elseif min(F(:)) >= 0 && max(F(:)) <= 255
    % Image in [0,255]
    Q_CI = std(F(:) / 255, 1);

else
    error('Input image must have values in either [0,1] or [0,255].');
end

end