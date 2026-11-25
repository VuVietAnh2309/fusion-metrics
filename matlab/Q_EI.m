function EI = Q_EI(img)
%Q_EI Calculate Edge Intensity (EI)
%
% Formula:
%
%       sqrt(sum_i sum_j (Sx(i,j)^2 + Sy(i,j)^2))
% EI = ------------------------------------------------
%                        M * N
%
% where:
%   Sx : horizontal Sobel gradient
%   Sy : vertical Sobel gradient
%   M  : number of rows
%   N  : number of columns
%
% Input:
%   img : grayscale or RGB image
%
% Output:
%   EI  : Edge Intensity value


    %% Convert RGB image to grayscale
    if ndims(img) == 3
        img = rgb2gray(img);
    end

    %% Convert to double
    img = double(img);

    %% Image size
    [M, N] = size(img);

    %% Sobel kernels
    hx = [ ...
         1  0 -1;
         2  0 -2;
         1  0 -1 ...
    ];

    hy = [ ...
         1  2  1;
         0  0  0;
        -1 -2 -1 ...
    ];

    %% Calculate gradients
    Sx = imfilter(img, hx, 'replicate', 'same', 'conv');
    Sy = imfilter(img, hy, 'replicate', 'same', 'conv');

    %% Calculate EI
    EI = sqrt(sum(Sx(:).^2 + Sy(:).^2)) / (M * N);

end