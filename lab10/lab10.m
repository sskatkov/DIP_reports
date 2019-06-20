%% Task 1. Inverse filter
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);
figure; imshow(I_corr); title('Motion & Noise');

I_corr = deconvwnr((I_corr), PSF);
figure; imshow(I_corr); title('Inverse filter (Motion & Noise)'); drawnow;

I_motion = deconvwnr((I_motion), PSF);
figure; imshow(I_motion); title('Inverse filter (Motion)'); drawnow;

%% Task 1. Wiener filter(1)
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);

noise  = double(I_corr) - double(I);
S_n = abs(fft2(noise)) .^ 2;

S_f = abs(fft2(I)) .^ 2;
NSR = sum (S_n(:)) / sum (S_f(:));

I_corr = deconvwnr((I_corr), PSF, NSR);
figure; imshow(I_corr); title('Wiener filter(1)'); drawnow;

%% Task 1. Wiener filter(2)
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);

[M,N] = size(I);

noise = double(I_corr) - double(I);
S_n = var(noise(:) - double(I(:)));
S_g = abs(fft2(I_corr) / sqrt(M*N)) .^ 2;
H = freqz2(PSF,M,N);
S = S_g - S_n;
S (S < 0) = 0;
S (S < 1) = 0;
S_f = S ./ (abs(H) .^ 2);
NSR = M * N * sum (S_n(:)) ./ sum (S_f(:));

I_corr = deconvwnr((I_corr), PSF, NSR);
figure; imshow(I_corr); title('Wiener filter(2)'); drawnow;

%% Task 1. Tikhonov regularization
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);
[M,N] = size(I);

I_corr = deconvreg(I_corr, PSF, 4.45 / 1000 * M * N, [1e-7 1e7]);
figure; imshow(I_corr); title('Tikhonov regularization');


%% Task 2. Inverse filter
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);
PSE = fspecial('gaussian',50, 10);
I_motion = edgetaper(I_motion, PSE);
I_corr = edgetaper(I_corr, PSE);

I_corr = deconvwnr((I_corr), PSF);
figure; imshow(I_corr); title('Inverse filtering (Motion & Noise)'); drawnow;

I_motion = deconvwnr((I_motion), PSF);
figure; imshow(I_motion); title('Inverse filtering (Motion)'); drawnow;

%% Task 2. Wiener filter(1)
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);
PSE = fspecial('gaussian',50, 10);
I_corr = edgetaper(I_corr, PSE);

noise  = double(I_corr) - double(I);
S_n = abs(fft2(noise)) .^ 2;

S_f = abs(fft2(I)) .^ 2;
NSR = sum (S_n(:)) / sum (S_f(:));

I_corr = deconvwnr((I_corr), PSF, NSR);
figure; imshow(I_corr); title('Wiener filter(1)'); drawnow;

%% Task 2. Wiener filter(2)
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);
PSE = fspecial('gaussian',50, 10);
I_corr = edgetaper(I_corr, PSE);

[M,N] = size(I);

noise = double(I_corr) - double(I);
S_n = var(noise(:) - double(I(:)));
S_g = abs(fft2(I_corr) / sqrt(M*N)) .^ 2;
H = freqz2(PSF,M,N);
S = S_g - S_n;
S (S < 0) = 0;
S (S < 1) = 0;
S_f = S ./ (abs(H) .^ 2);
NSR = M * N * sum (S_n(:)) ./ sum (S_f(:));

I_corr = deconvwnr((I_corr), PSF, NSR);
figure; imshow(I_corr); title('Wiener filter(2)'); drawnow;

%% Task 2. Tikhonov regularization
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);
PSE = fspecial('gaussian',50, 10);
I_corr = edgetaper(I_corr, PSE);
[M,N] = size(I);

I_corr = deconvreg(I_corr, PSF, 4.45 / 1000 * M * N, [1e-7 1e7]);
figure; imshow(I_corr); title('Tikhonov regularization'); drawnow;

%% Task 3
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
PSE = fspecial('gaussian', 50, 10);
I_motion = edgetaper(I_motion, PSE);
figure; imshow(I_motion); title('Motion'); drawnow;

G = deconvlucy(I_motion, PSF, 50);
figure; imshow(G); title('Richardson-Lucy deconvolution'); drawnow;

%% Task 4
clc; close all; clear;
I = imread('up.png'); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
figure; imshow(I_motion); title('Motion'); drawnow;

G = deconvlucy(I_motion, PSF, 50);
figure; imshow(G); title('Richardson–Lucy deconvolution'); drawnow;

%% Task 5
clc; close all; clear;
I = im2double(imread('up.png')); figure; imshow(I); title('Original'); drawnow;
PSF = fspecial('motion', 10, 45);
I_motion = imfilter(I, PSF, 'replicate');
I_corr = imnoise(I_motion, 'gaussian', 0, 0.005);
PSE = fspecial('gaussian', 50, 10);
figure; imshow(I_corr); title('Motion & Noise'); drawnow;

G = deconvlucy(I_corr, PSF, 50, 0.12);
figure; imshow(G); title('Richardson-Lucy deconvolution'); drawnow;

%%
% 
%  Все методы кроме инверсной фильтрации показали нормальные результаты. Лучших результатов позволяет добиться метод Люси-Ричардсона, но и он при возникновении шума на изображении дает далеко не идеальный результат,
%  кроме того, при восстановлении практически необходимо использовать edgetraper -- в таком случае не возникает "звон" по краям изображения. При восстановлении необходимо эмпирически подбирать параметры метода.
% 
