% Q_TE -- Cvejic, Canagarajah & Bull's Tsallis-entropy MI fusion metric
% (N. Cvejic et al., "Image fusion metric based on mutual information
% and Tsallis entropy," Electronics Letters 42(11), 2006, eq. 9-10 --
% see Paper_Metrics/cvejic2006_Q_TE.pdf, also cited by Liu et al. 2012 Sec.
% 2.1.2 as this metric's origin).
%
% Corrected 2026-08-28: the formula here used to match Liu et al. 2012's
% own eq. 5 verbatim -- but that eq. transcribes Cvejic's denominator
% wrong. Cvejic's actual eq. 9 is (p_F(f)*p_X(x))^(alpha-1) -- BOTH
% marginals raised to the SAME (q-1) power (symmetric in X and Y, as any
% proper divergence-from-independence measure should be, since MI is
% symmetric). Liu 2012's eq. 5 instead writes h_F(j)*h_A(i)^(q-1) --
% only ONE marginal exponentiated -- an asymmetric formula that does not
% match its own cited source. Re-derived directly from Cvejic's eq. 9
% after the user supplied the original paper. Default q also corrected
% to 1.85: Cvejic's paper explicitly states "extensive subjective
% experiments... alpha=1.85 was chosen" as the empirically optimal
% value -- q=1.5 (the old default) has no cited justification anywhere.
function Q = Q_TE(A, B, F, q)

if nargin < 4
    q = 1.85; % default (Cvejic 2006's own empirically-tuned value)
end

grey_level = 256;

A = uint8(A);
B = uint8(B);
F = uint8(F);

I_AF = Tsallis_MI(A, F, q, grey_level);
I_BF = Tsallis_MI(B, F, q, grey_level);

Q = I_AF + I_BF;

end

function Iq = Tsallis_MI(X, Y, q, grey_level)

X = double(X);
Y = double(Y);

% ===== Joint histogram =====
joint_hist = accumarray([X(:)+1, Y(:)+1], 1, [grey_level grey_level]);

p_xy = joint_hist / sum(joint_hist(:));

% ===== Marginal =====
p_x = sum(p_xy, 2);   % h_A(i)
p_y = sum(p_xy, 1);   % h_F(j)

[i,j] = find(p_xy > 0);

p = p_xy(sub2ind(size(p_xy), i, j));
px = p_x(i);
py = p_y(j)';

% ===== Tsallis MI (Cvejic 2006 eq. 9 -- symmetric marginals) =====
Iq = (1/(1-q)) * (1 - sum( (p.^q) ./ ((px.*py).^(q-1)) ));

end