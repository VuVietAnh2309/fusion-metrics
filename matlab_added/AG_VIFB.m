function res = AG_VIFB(F)
% Average gradient, VIFB definition (VIFB toolbox, metricsAvg_gradient.m;
% X. Zhang, P. Ye, G. Xiao, "VIFB: A visible and infrared image fusion benchmark", CVPRW 2020).
% NOT the reference Q_AG.m (which uses forward differences on the [0,1] scale).
% Input: grayscale fused image, values 0..255. Output on the 0..255 scale.
F = double(F);
[r, c] = size(F);
[dzdx, dzdy] = gradient(F);
s = sqrt((dzdx .^ 2 + dzdy .^ 2) ./ 2);
res = sum(sum(s)) / ((r - 1) * (c - 1));
end
