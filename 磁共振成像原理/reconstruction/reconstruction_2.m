%% a
accelerationFactor = 4;
n = 256;

k_len = n;

acceleration_factor = 4;
numSelected = k_len/acceleration_factor;
distances = abs((1:k_len) - k_len/2);
weights = 1./(distances + 1);
selectedPoints = datasample(1:k_len, numSelected, 'Weights', weights, 'Replace', false);
amplitudes = zeros(1, k_len);
amplitudes(selectedPoints) = 1;
matrixAmplitudes = repmat(amplitudes, k_len, 1);

UR4 = zeros(256,256);
imgSize = 256;

prm = randperm(size(M,1));
density = rand(imgSize, 1);
UR4(:,prm(:,1:size(M,1)/accelerationFactor)) = 1;
UVDR4 = matrixAmplitudes; 

PSF_UR4 = fftshift(fft2(ifftshift(UR4)));
PSF_UVDR4 = fftshift(fft2(ifftshift(UVDR4)));

aliased_UR4 = ifft2(UR4 .* fft2(M));
aliased_UVDR4 = ifft2(UVDR4 .*fft2(M));

figure;
subplot(2, 2, 1);
imshow(UR4, []);
title('Random Undersampling (UR4)');
subplot(2, 2, 2);
imshow(UVDR4, []);
title('Variable Density Random Undersampling (UVDR4)');
subplot(2, 2, 3);
imshow(abs(PSF_UR4), []);
title('PSF for UR4');
subplot(2, 2, 4);
imshow(abs(PSF_UVDR4), []);
title('PSF for UVDR4');

figure;
subplot(1,2,1);
imshow(abs(aliased_UR4),[]);
title('UR4')
subplot(1,2,2);
imshow(abs(aliased_UVDR4),[]);
title('UVDR4')

%% b
% L1-magic
n = 512; % 信号长度
k = 50; % 稀疏度
x = zeros(n, 1);
x(randperm(n, k)) = randn(k, 1);

% 生成观测测量 b（添加噪声模拟测量）
A_UR4 = randn(100, n); % 观测矩阵
sigma = 0.1; % 噪声水平
b = A_UR4 * x + sigma * randn(100, 1);

% 使用 L1-Magic 进行基础的 Basis Pursuit DeNoising (BPDN) 重建
% 这里使用 L1-magic 中的 l1eq_pd 函数

% 设置参数
lambda = 0.1; % 正则化参数

% 调用 l1eq_pd 函数进行重建
reconstructed_x = l1eq_pd(x, A_UR4, [], b, lambda);

% 显示结果
figure;
subplot(2, 1, 1);
stem(x, 'r', 'filled');
title('True Sparse Signal');

subplot(2, 1, 2);
stem(reconstructed_x, 'b', 'filled');
title('Reconstructed Signal using L1-Magic');

%% c
f = M;
mask_UR4 = UR4;
A_UR4 = @(x)  masked_FFT(x,mask_UR4);
AT_UR4 = @(x) masked_FFT_t(x,mask_UR4);
y_UR4 = A_UR4(f);

mask_UVDR4 = UVDR4;
A_UVDR4 = @(x)  masked_FFT(x,mask_UVDR4);
AT_UVDR4 = @(x) masked_FFT_t(x,mask_UVDR4);
y_UVDR4 = A_UVDR4(f);
% denoising function;
tv_iters = 10;
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);
Phi = @(x) TVnorm(x);

tau = 0.001;
[x_twist_UR4,~,objective_UR4, times_UR4,~,mse_UR4] = TwIST(y_UR4,A_UR4,tau,'Lambda', 1e-4,'AT', AT_UR4, 'Psi', Psi, 'Phi',Phi, 'True_x', f, 'MaxiterA', 1000, 'StopCriterion',1, 'ToleranceA',1e-4, 'Verbose', 1);
[x_twist_UVDR4,~,objective_UVDR4, times_UVDR4,~,mse_UVDR4] = TwIST(y_UVDR4,A_UVDR4,tau,'Lambda', 1e-4,'AT', AT_UVDR4, 'Psi', Psi, 'Phi',Phi, 'True_x', f, 'MaxiterA', 1000, 'StopCriterion',1, 'ToleranceA',1e-4, 'Verbose', 1);
figure;
subplot(1,3,1);
imshow(abs(f));
colormap(gray);
axis image;
axis off;
title('Original');
subplot(1,3,2);
imagesc(x_twist_UR4);
colormap(gray);
axis image;
axis off;
title('Estimate_{UR4}');
subplot(1,3,3);
imagesc(x_twist_UVDR4);
axis image;
axis off;
title('Estimate_{UVDR4}');
%% d




