function out = blkproc(A, blk, border, fun)
% Stand-in used ONLY when running Q_CV.m under GNU Octave (the image package no longer ships
% blkproc). MATLAB's blkproc(A, [m n], [0 0], fun) zero-pads A at the bottom/right to a multiple
% of the block size, applies fun to every m-by-n block and returns the matrix of results
% (fun returns a scalar here).
if any(border ~= 0)
    error('shim blkproc: only a [0 0] border is provided');
end
m = blk(1); n = blk(2);
[H, W] = size(A);
Hp = ceil(H / m) * m; Wp = ceil(W / n) * n;
P = zeros(Hp, Wp, class(A));
P(1:H, 1:W) = A;
out = zeros(Hp / m, Wp / n);
for i = 1:Hp / m
    for j = 1:Wp / n
        out(i, j) = fun(P((i-1)*m+1:i*m, (j-1)*n+1:j*n));
    end
end
end
