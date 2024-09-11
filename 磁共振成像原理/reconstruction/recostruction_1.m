%% a
Nc = 8;
mc = zeros(size(M, 1), size(M, 2), Nc);

for i = 1:Nc
    mc(:,:,i) = C(:,:,i) .* M;
end
figure;
for i = 1:Nc
   
    subplot(1, Nc, i);
    imshow(abs(mc(:,:,i)), []);
    title(['Coil Image ' num2str(i)]);
end
 
rss_combined = sqrt(sum(abs(mc).^2, 3));

weighted_combined = sum(mc .* conj(C), 3);

figure;


subplot(1, 2, 1);
imshow(abs(rss_combined), []);
title('Root-Sum-of-Square (RSS) Combined Image');

subplot(1, 2, 2);
imshow(abs(weighted_combined), []);
title('Weighted Coil Sensitivity Combined Image');

%% b
imgSize = 256; 

accFactors = [3, 7];

U3 = zeros(imgSize, imgSize);
U7 = zeros(imgSize, imgSize);
UR3 = zeros(imgSize, imgSize);
UR7 = zeros(imgSize, imgSize);

prm = randperm(size(M,1));
U3(:,1:3:end) = 1;
U7(:,1:7:end) = 1;
UR3(:,prm(:,1:size(M,1)/3)) = 1;
UR7(:,prm(:,1:size(M,1)/7)) = 1;

PSF_U3 = fftshift(fft2(ifftshift(U3)));
PSF_U7 = fftshift(fft2(ifftshift(U7)));
PSF_UR3 = fftshift(fft2(ifftshift(UR3)));
PSF_UR7 = fftshift(fft2(ifftshift(UR7)));

figure;
subplot(2, 4, 1), imshow(U3, []), title('U3 Pattern');
subplot(2, 4, 2), imshow(U7, []), title('U7 Pattern');
subplot(2, 4, 3), imshow(UR3, []), title('UR3 Pattern');
subplot(2, 4, 4), imshow(UR7, []), title('UR7 Pattern');
subplot(2, 4, 5), imshow(log(abs(PSF_U3)+1), []), title('PSF of U3');
subplot(2, 4, 6), imshow(log(abs(PSF_U7)+1), []), title('PSF of U7');
subplot(2, 4, 7), imshow(log(abs(PSF_UR3)+1), []), title('PSF of UR3');
subplot(2, 4, 8), imshow(log(abs(PSF_UR7)+1), []), title('PSF of UR7');
%% c
Nc = size(C, 3);

aliased_U3 = zeros(size(M, 1), size(M, 2), Nc);
aliased_U7 = zeros(size(M, 1), size(M, 2), Nc);
aliased_UR3 = zeros(size(M, 1), size(M, 2), Nc);
aliased_UR7 = zeros(size(M, 1), size(M, 2), Nc);

for i = 1:Nc
    Ci_m = C(:, :, i) .* M;

    F_Ci_m = fft2(Ci_m);

    aliased_U3(:, :, i) = ifft2(U3 .* F_Ci_m);
    aliased_U7(:, :, i) = ifft2(U7 .* F_Ci_m);
    aliased_UR3(:, :, i) = ifft2(UR3 .* F_Ci_m);
    aliased_UR7(:, :, i) = ifft2(UR7 .* F_Ci_m);
end

for i = 1:Nc
    figure;
    subplot(2, 2, 1), imshow(abs(aliased_U3(:, :, i)), []), title(['Coil ', num2str(i), ' with U3']);
    subplot(2, 2, 2), imshow(abs(aliased_U7(:, :, i)), []), title(['Coil ', num2str(i), ' with U7']);
    subplot(2, 2, 3), imshow(abs(aliased_UR3(:, :, i)), []), title(['Coil ', num2str(i), ' with UR3']);
    subplot(2, 2, 4), imshow(abs(aliased_UR7(:, :, i)), []), title(['Coil ', num2str(i), ' with UR7']);
end

%% d
E3 = Cartesian_SENSE(U3,C);
b3 = E3 * M;
rho3 = gradientDescent(b3,E3,1000);

E7 = Cartesian_SENSE(U7,C);
b7 = E7 * M;
rho7 = gradientDescent(b7,E7,1000);

ER3 = Cartesian_SENSE(UR3,C);
bR3 = ER3 * M;
rhoR3 = gradientDescent(bR3,ER3,1000);

ER7 = Cartesian_SENSE(UR7,C);
bR7 = ER7 * M;
rhoR7 = gradientDescent(bR7,ER7,1000);

figure;
subplot(2,2,1)
imshow(abs(rho3(:,:,end)),[]);title('Gradient_{U_3}');
subplot(2,2,2)
imshow(abs(rho7(:,:,end)),[]);title('Gradient_{U_7}');
subplot(2,2,3)
imshow(abs(rhoR3(:,:,end)),[]);title('Gradient_{U_R3}');
subplot(2,2,4)
imshow(abs(rhoR7(:,:,end)),[]);title('Gradient_{U_R7}');
%% e
iteration = 1:10:100;
sigma_U3 = zeros(1, length(iteration));
sigma_U7 = zeros(1, length(iteration));
sigma_UR3 = zeros(1, length(iteration));
sigma_UR7 = zeros(1, length(iteration));

for j = 1:length(iteration)
    i = iteration(j);

    rho3 = abs(gradientDescentMSE(b3, E3, i));
    rho7 = abs(gradientDescentMSE(b7, E7, i));
    rhoR3 = abs(gradientDescentMSE(bR3, ER3, i));
    rhoR7 = abs(gradientDescentMSE(bR7, ER7, i));

    error_U3 = zeros(256, 256);
    error_U7 = zeros(256, 256);
    error_UR3 = zeros(256, 256);
    error_UR7 = zeros(256, 256);

    for m = 1:256
        for n = 1:256
            error_U3(m, n) = (rho3(m, n) - M(m, n))^2;
            error_U7(m, n) = (rho7(m, n) - M(m, n))^2;
            error_UR3(m, n) = (rhoR3(m, n) - M(m, n))^2;
            error_UR7(m, n) = (rhoR7(m, n) - M(m, n))^2;
        end
    end

    sigma_U3(j) = sum(error_U3(:)) / (256 * 256);
    sigma_U7(j) = sum(error_U7(:)) / (256 * 256);
    sigma_UR3(j) = sum(error_UR3(:)) / (256 * 256);
    sigma_UR7(j) = sum(error_UR7(:)) / (256 * 256);
    
end



figure;
plot(iteration, sigma_U3, '-o', 'DisplayName', 'Uniform U3');
hold on;
plot(iteration, sigma_U7, '-s', 'DisplayName', 'Uniform U7');
hold on;
plot(iteration, sigma_UR3, '-^', 'DisplayName', 'Random UR3');
hold on;
plot(iteration, sigma_UR7, '-d', 'DisplayName', 'Random UR7');
hold off;

xlabel('Iteration');
ylabel('Sigma');
title('Relationship between Sigma and Iteration');
legend('show');
grid on;
hold off;

