function [f1, f2] = freqspace(n, flag)
% Stand-in used ONLY when running Q_CV.m / Q_CB.m under GNU Octave, which does not implement
% freqspace. MATLAB definition: for an m-by-n matrix (n = [m n])
%   f1 = ((0:n(2)-1) - floor(n(2)/2)) * (2/n(2))     % columns
%   f2 = ((0:n(1)-1) - floor(n(1)/2)) * (2/n(1))     % rows
% i.e. (-N/2:N/2-1)*2/N for even N and (-(N-1)/2:(N-1)/2)*2/N for odd N, as documented;
% 'meshgrid' returns meshgrid(f1, f2).
if numel(n) == 1
    n = [n n];
end
f1 = ((0:n(2)-1) - floor(n(2)/2)) * (2/n(2));
f2 = ((0:n(1)-1) - floor(n(1)/2)) * (2/n(1));
if nargin > 1 && strcmpi(flag, 'meshgrid')
    [f1, f2] = meshgrid(f1, f2);
end
end
