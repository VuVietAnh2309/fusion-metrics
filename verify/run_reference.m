% Compute the reference values of the self-test with the ORIGINAL .m files (GNU Octave).
%   octave-cli run_reference.m     (reads reference_cases.mat, writes reference_values.mat)
% Shims (octave_shims/): dwt2 (Haar, absent in Octave) and imgradientxy (MATLAB boundary
% convention). Q_FMI is evaluated with feature 'none' on stored Sobel edge maps, because
% Octave's edge() thins differently from MATLAB's.
warning('off', 'all');
pkg load image
pkg load signal
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'matlab'));
addpath(fullfile(here, '..', 'matlab_added'));
addpath(fullfile(here, 'octave_shims'));
S = load(fullfile(here, 'reference_cases.mat'));
names = {'Q_MLI','Q_AG','Q_En','Q_SF','Q_PSNR','Q_CC','Q_SCD','Q_MI','Q_ABF','Q_VIFF', ...
         'Q_M','Q_NCIE','Q_TE_Nava','Q_FSIM','AG_VIFB','VIF_Sheikh','Q_FMI_core'};
K = numel(S.names);
R = nan(K, numel(names));
for k = 1:K
    A = S.(sprintf('A%d', k)); B = S.(sprintf('B%d', k)); F = S.(sprintf('F%d', k));
    R(k,1) = Q_MLI(F);  R(k,2) = Q_AG(F);  R(k,3) = Q_En(F);
    R(k,4) = Q_SF(A,B,F);  R(k,5) = Q_PSNR(A,B,F);  R(k,6) = Q_CC(A,B,F);
    R(k,7) = Q_SCD(A,B,F); R(k,8) = Q_MI(A,B,F,256); R(k,9) = Q_ABF(A,B,F);
    R(k,10) = Q_VIFF(A,B,F); R(k,11) = Q_M(A,B,F); R(k,12) = Q_NCIE(A,B,F);
    R(k,13) = Q_TE_Nava(A,B,F); R(k,14) = Q_FSIM(A,B,F);
    R(k,15) = AG_VIFB(F); R(k,16) = VIF_Sheikh(A,B,F);
    R(k,17) = Q_FMI(double(S.(sprintf('EA%d', k))), double(S.(sprintf('EB%d', k))), ...
                    double(S.(sprintf('EF%d', k))), 'none', 3);
    printf('case %d/%d done\n', k, K); fflush(stdout);
end
save('-v7', fullfile(here, 'reference_values.mat'), 'R', 'names');
