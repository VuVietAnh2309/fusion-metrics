function [cA, cH, cV, cD] = dwt2(x, wname)
% Minimal stand-in for MATLAB's dwt2(x, 'haar') (Wavelet Toolbox), used ONLY to run the
% reference Q_M.m under Octave, which has no dwt2. Default MATLAB extension mode 'sym'.
% For the Haar filter this equals non-overlapping 2x2 sums/differences scaled by 1/2.
if ~strcmp(wname, 'haar')
    error('shim dwt2: only haar is supported');
end
x = double(x);
if mod(size(x,1),2), x = [x; x(end,:)]; end
if mod(size(x,2),2), x = [x, x(:,end)]; end
s = 1/sqrt(2);
lo = (x(:,1:2:end) + x(:,2:2:end)) * s;
hi = (x(:,1:2:end) - x(:,2:2:end)) * s;
cA = (lo(1:2:end,:) + lo(2:2:end,:)) * s;
cH = (lo(1:2:end,:) - lo(2:2:end,:)) * s;
cV = (hi(1:2:end,:) + hi(2:2:end,:)) * s;
cD = (hi(1:2:end,:) - hi(2:2:end,:)) * s;
end
