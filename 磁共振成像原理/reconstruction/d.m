Nc = 8;
U3 = zeros(imgSize, imgSize);
U7 = zeros(imgSize, imgSize);
UR3 = zeros(imgSize, imgSize);
UR7 = zeros(imgSize, imgSize);

% Generate uniform undersampling patterns
for i = 1:imgSize
    if mod(i, accFactors(1)) == 0
        U3(i, :) = 1;
    end
    if mod(i, accFactors(2)) == 0
        U7(i, :) = 1;
    end
end

% Generate random undersampling patterns
UR3(randperm(numel(UR3), round(imgSize^2/accFactors(1)))) = 1;
UR7(randperm(numel(UR7), round(imgSize^2/accFactors(2)))) = 1;

for i = 1:Nc
    b(i) = ifft(U3*fft(C(:,:,i) .* M));
    imshow(abs(bi));
end

