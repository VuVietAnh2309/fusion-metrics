% Q_En -- Shannon Entropy. NOT part of Liu et al. 2012's 12-metric
% taxonomy (single-image, no-reference statistic).
% Copied unchanged from Metrics_Matlab/Q_En.m.
function E = Q_En(F)
% Q_CI - Contrast Index

F = double(F);

if min(F(:)) >= 0 && max(F(:)) <= 1
    % Image in [0,1]
    E = entropy(F);

elseif min(F(:)) >= 0 && max(F(:)) <= 255
    % Image in [0,255]
    E = entropy(F./255);

else
    error('Input image must have values in either [0,1] or [0,255].');
end

end