function res = Q_NMI(im1, im2, fim)
% Q_NMI -- Liu et al. 2012's "Normalized Mutual Information" fusion metric
% (Section 2.1.1, eq. 4), using Hossny et al.'s normalized definition:
%
%   Q_NMI = 2*[ MI(A,F)/(H(A)+H(F)) + MI(B,F)/(H(B)+H(F)) ]
%
% Was a thin wrapper calling ../../_deps_matlab/Q_NMI.m (the shared,
% unchanged multi-variant dependency, sw=1..3) with sw=1 fixed -- but
% after the 2026-08-27 rename (Q_MI -> Q_NMI, to free up "Q_MI" for the
% Qu 2006 metric elsewhere in this group), THIS file is now ALSO named
% Q_NMI.m, which would shadow-collide with ../../_deps_matlab/Q_NMI.m on
% the MATLAB path (whichever is added last wins, and a bare `Q_NMI(...)`
% call from inside this very file would then recurse into itself instead
% of reaching the dependency). Fixed by inlining the sw=1 branch's actual
% math directly (still calls the shared, still-unchanged normalize1.m/
% mutual_info.m helpers in ../../_deps_matlab -- those names don't
% collide with anything here) instead of calling Q_NMI(...,1) by name.
% Requires normalize1.m and mutual_info.m on the MATLAB path (see
% ../../_deps_matlab).
%
% NOTE: this is NOT the same as Q_MI.m in this same folder (Qu 2006's
% unnormalized MI, a different metric not covered by this paper) --
% different formula, similar author-chosen name. Don't confuse the two.

im1 = normalize1(im1);
im2 = normalize1(im2);
fim = normalize1(fim);

[I_fx, ~, H_x, H_f1] = mutual_info(im1, fim);
[I_fy, ~, H_y, H_f2] = mutual_info(im2, fim);

res = 2 * (I_fx / (H_f1 + H_x) + I_fy / (H_f2 + H_y));
end


%---------------------------------------------------------------------------
% Inlined from Metrics_Full/_deps_matlab/mutual_info.m, normalize1.m so this file has no
% sibling-file dependency (mutual_info, normalize1).
%---------------------------------------------------------------------------

function [MI,H_xy,H_x,H_y]=mutual_info(im1,im2)

% function [MI,H_xy,H_x,H_y]=mutual_info(im1,im2)
% 
% This function is caculate the mutual information of two input images.
% im1   -- input image one;
% im2   -- input image two;
%
% MI    -- mutual information;
%
%
% Note: The input images need to be in the range of 0-255. (see function:
% normalize1.m)
%
% Z. Liu @ NRCC [July 17, 2009]

im1=double(im1);
im2=double(im2);

[hang,lie]=size(im1);
count=hang*lie;
N=256;

%% caculate the joint histogram
h=zeros(N,N);

for i=1:hang
    for j=1:lie
        % in this case im1->x (row), im2->y (column)
        h(im1(i,j)+1,im2(i,j)+1)=h(im1(i,j)+1,im2(i,j)+1)+1;
    end
end

%% marginal histogram

% this operation converts histogram to probability
%h=h./count;
h=h./sum(h(:));

im1_marg=sum(h);    % sum each column for im1
im2_marg=sum(h');   % sum each row for im2

H_x=0;   % for im1
H_y=0;   % for im2

%for i=1:count
%for i=1:N
%    if (im1_marg(i)>eps)
%        % entropy for image1
%        H_x=H_x+(-im1_marg(i)*(log2(im1_marg(i))));
%    end
%    if (im2_marg(i)>eps)
%        % entropy for image2
%        H_y=H_y+(-im2_marg(i)*(log2(im2_marg(i))));
%    end
%end

H_x=-sum(im1_marg.*log2(im1_marg+(im1_marg==0)));
H_y=-sum(im2_marg.*log2(im2_marg+(im2_marg==0)));


% joint entropy
H_xy=-sum(sum(h.*log2(h+(h==0))));

% mutual information
MI=H_x+H_y-H_xy;
end

function RES=normalize1(data)

% function RES=normalize1(data)
%
% This function is to NORMALIZE the data. 
% The data will be in the interval 0-255 (gray level) and pixel value has
% been rounded to an integer.
% 
% See also: normalize.m 
%
% Z. Liu @NRCC (Aug 24, 2009)

data=double(data);
da=max(data(:));
xiao=min(data(:));
if (da==0 & xiao==0)
    RES=data;
else
    newdata=(data-xiao)/(da-xiao);
    RES=round(newdata*255);
end
end
