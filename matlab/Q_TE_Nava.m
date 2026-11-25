% Q_TE_Nava -- Nava, Cristobal & Escalante-Ramirez's NORMALIZED Tsallis
% mutual-information fusion metric (R. Nava, G. Cristobal, B. Escalante-
% Ramirez, "Non-reference image fusion evaluation procedure based on
% mutual information and a generalized entropy measure," Proc. SPIE
% 6592, 659215, 2007 -- see Paper_Metrics/nava2007_Q_TE.pdf, eq. 11-16).
% Cited by Liu et al. 2012 alongside Cvejic 2006 for Q_TE (Sec. 2.1.2);
% Nava's own contribution is the NORMALIZED version -- Liu 2012's eq. 7
% IS this paper's eq. 16.
%
% Nava's eq. 11-12 for the (un-normalized) Tsallis MI between two images
% matches Cvejic 2006's eq. 9 exactly (both marginals sharing the same
% (q-1) exponent, symmetric) -- an independent confirmation that the
% symmetric form is correct, and that Liu 2012's eq. 5 (asymmetric) is a
% transcription error (see Q_TE.m's comment for the full story).
%
% Nava's actual proposed metric normalizes the raw sum
% I_q(F,A)+I_q(F,B) (Cvejic's M^q_FAB, eq. 13) by the Tsallis JOINT
% entropy of the two SOURCE images, H_q(A,B) = H_q(A)+H_q(B)-I_q(A,B)
% (eq. 15) -- "how much of the total information available across both
% sources made it into F", bounded and comparable across image pairs,
% unlike Cvejic's raw (unbounded) sum:
%
%   NM_q(F,A,B) = [I_q(F,A) + I_q(F,B)] / [H_q(A) + H_q(B) - I_q(A,B)]
%
% q has no single value Nava 2007 itself calls "the" default -- their
% own experiments use three different values (0.2624, 0.43137, 0.6028),
% all < 1. Defaulting to Cvejic 2006's q=1.85 here (as this project's
% sibling Q_TE.m does) was tried and rejected: verified empirically that
% the denominator H_q(A)+H_q(B)-I_q(A,B) collapses toward 0 as q grows
% past 1 (Tsallis entropy itself shrinks fast for q>1 on a 256-level
% real-image histogram), pushing the ratio past 10 -- far outside the
% ~0-1 range every result in Nava's own paper reports. Default set to
% q=0.43137 -- the value subsequent image-fusion literature settled on
% when reusing Nava's normalized metric (this project's user,
% 2026-08-28). Pass an explicit q to reproduce one of Nava's other two
% reported experiments (0.2624, 0.6028) instead.
%
% Written from scratch (no reference MATLAB toolbox implementation of
% this normalized variant was found anywhere in this project) directly
% from Nava 2007's equations. Cross-verified against Q_TE_Nava.py (an
% independent Python implementation of the same equations) rather than
% against a pre-existing third-party reference, since none exists for
% this specific normalized metric.
function NM = Q_TE_Nava(A, B, F, q)

if nargin < 4
    q = 0.43137; % default (value settled on by subsequent image-fusion literature)
end

grey_level = 256;

A = uint8(A);
B = uint8(B);
F = uint8(F);

I_FA = Tsallis_MI(F, A, q, grey_level);
I_FB = Tsallis_MI(F, B, q, grey_level);
H_A  = Tsallis_Entropy(A, q, grey_level);
H_B  = Tsallis_Entropy(B, q, grey_level);
I_AB = Tsallis_MI(A, B, q, grey_level);

NM = (I_FA + I_FB) / (H_A + H_B - I_AB);

end


function H = Tsallis_Entropy(X, q, grey_level)
% H_q(X), eq. 9 of Nava 2007.
X = double(X);
counts = accumarray(X(:) + 1, 1, [grey_level 1]);
p = counts / sum(counts);
p = p(p > 0);
H = (1 / (1 - q)) * (sum(p .^ q) - 1);
end


function Iq = Tsallis_MI(X, Y, q, grey_level)
% I_q(X;Y), eq. 11-12 of Nava 2007 -- matches Cvejic 2006 eq. 9
% (symmetric marginals).
X = double(X);
Y = double(Y);

joint_hist = accumarray([X(:)+1, Y(:)+1], 1, [grey_level grey_level]);
p_xy = joint_hist / sum(joint_hist(:));

p_x = sum(p_xy, 2);
p_y = sum(p_xy, 1);

[i, j] = find(p_xy > 0);
p = p_xy(sub2ind(size(p_xy), i, j));
px = p_x(i);
py = p_y(j)';

Iq = (1 / (1 - q)) * (1 - sum((p .^ q) ./ ((px .* py) .^ (q - 1))));
end
