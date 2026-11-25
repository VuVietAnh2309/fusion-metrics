function res=Q_CV(im1,im2,fim)
% Q_CV -- Liu et al. 2012's "Chen-Varshney Metric" (Sec. 2.4.1, eq. 47-48).
%
% Ref: H. Chen and P.K. Varshney, "A human perception inspired quality
% metric for image fusion based on regional information," Information
% Fusion, vol. 8, pp. 193-207, 2007.
%
% Faithful port of metricChen.m (Z. Liu @NRCC, July 2009). Renamed to
% match this project's Q_<NAME>(A,B,F) convention. Uses the
% Mannos-Sakrison CSF filter (the paper's default choice) via the
% author's own frequency-domain filtering; Daly and Ahumada filter
% variants are computed but unused in the original, kept here in case
% they're wanted later. Needs normalize1.m (see ../../_deps_matlab).

alpha_c=1;
alpha_s=0.685;
f_c=97.3227;
f_s=12.1653;

% local window size 16 x 16
windowSize=16;

% alpha = 1, 2, 3, 4, 5, 10, 15. This value is adjustable.
alpha=5;

im1=double(im1);
im2=double(im2);
fim=double(fim);

im1=normalize1(im1);
im2=normalize1(im2);
fim=normalize1(fim);

%% Step 1: extract edge information

flt1=[-1 0 1 ; -2 0 2 ; -1 0 1];
flt2=[-1 -2 -1; 0 0 0; 1 2 1];

fuseX=filter2(flt1,fim,'same'); %#ok<NASGU>
fuseY=filter2(flt2,fim,'same'); %#ok<NASGU>

img1X=filter2(flt1,im1,'same');
img1Y=filter2(flt2,im1,'same');
im1G=sqrt(img1X.*img1X+img1Y.*img1Y);

img2X=filter2(flt1,im2,'same');
img2Y=filter2(flt2,im2,'same');
im2G=sqrt(img2X.*img2X+img2Y.*img2Y);

%% calculate the local region saliency

[hang,lie]=size(im1);

fun=@(x) sum(sum(x.^alpha));

ramda1=blkproc(im1G, [windowSize windowSize],[0 0],fun);
ramda2=blkproc(im2G, [windowSize windowSize],[0 0],fun);

%% similarity measurement

f1=im1-fim;
f2=im2-fim;

% filtering with CSF filters
% first build the three filters in frequency domain:

[u,v]=freqspace([hang,lie],'meshgrid');

% length 3
u=lie/8*u; v=hang/8*v;

r=sqrt(u.^2+v.^2);

% Mannos-Skarison's filter
theta_m=2.6*(0.0192+0.144*r).*exp(-(0.144*r).^1.1);

% second, filter the image in frequency domain
ff1=fft2(f1);
ff2=fft2(f2);

% here filter 1 (Mannos-Sakrison) is used, matching the paper's default.
Df1=ifft2(ifftshift(fftshift(ff1).*theta_m));
Df2=ifft2(ifftshift(fftshift(ff2).*theta_m));

fun2=@(x) mean2(x.^2);
D1=blkproc(Df1, [windowSize windowSize],[0 0],fun2);
D2=blkproc(Df2, [windowSize windowSize],[0 0],fun2);

%% global quality

Q=sum(sum(ramda1.*D1+ramda2.*D2))/sum(sum(ramda1+ramda2));

res=Q;

end


%---------------------------------------------------------------------------
% Inlined from Metrics_Full/_deps_matlab/normalize1.m so this file has no
% sibling-file dependency (normalize1).
%---------------------------------------------------------------------------

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
