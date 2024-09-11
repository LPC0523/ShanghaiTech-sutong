% 图像大小
n = 256;

% 生成高斯分布的权重矩阵
[X, Y] = meshgrid(1:n, 1:n);
center = [n/2, n/2]; % 高斯分布的中心位置
sigma = 0.1 * n; % 控制高斯分布的标准差，可以根据需求调整

% 计算高斯分布的权重
weights = exp(-((X - center(1)).^2 + (Y - center(2)).^2) / (2 * sigma^2));

% 归一化权重，确保概率之和为1
weights = weights / sum(weights(:));

% 随机生成图像（模拟 k 空间中的数据）
imageInKSpace = fft2(randn(n));

% 将权重应用在 k 空间
weightedImageInKSpace = imageInKSpace .* fftshift(weights);

% 反变换得到带有权重的图像
weightedImageInSpace = ifft2(weightedImageInKSpace);

% 显示结果
figure;

subplot(1, 3, 1);
imagesc(abs(imageInKSpace));
axis equal;
title('Original k-Space');

subplot(1, 3, 2);
imagesc(abs(weightedImageInKSpace));
axis equal;
title('Weighted k-Space');

subplot(1, 3, 3);
imagesc(abs(weightedImageInSpace));
axis equal;
title('Weighted Image in Space');
