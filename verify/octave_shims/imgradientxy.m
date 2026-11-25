function [Gx, Gy] = imgradientxy(I, method)
% Stand-in used ONLY when running the reference Q_AG.m under GNU Octave.
% Octave's image package computes 'intermediate' differences with conv2(...,'same'), which gives
% -I at the last column/row; MATLAB's imgradientxy sets the last column/row to 0:
%   Gx(:,1:end-1) = diff(I,1,2);  Gy(1:end-1,:) = diff(I,1,1).
% Values computed in MATLAB confirm this convention.
if nargin < 2 || ~any(strcmpi(method, {'intermediate', 'intermediatedifference'}))
    error('shim imgradientxy: only the ''intermediate'' method is provided');
end
I = double(I);
Gx = zeros(size(I));
Gy = zeros(size(I));
if size(I, 2) > 1
    Gx(:, 1:end-1) = diff(I, 1, 2);
end
if size(I, 1) > 1
    Gy(1:end-1, :) = diff(I, 1, 1);
end
end
