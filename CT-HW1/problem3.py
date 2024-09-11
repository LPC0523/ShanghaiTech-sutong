import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage.metrics import mean_squared_error, peak_signal_noise_ratio, structural_similarity

noisy_image = cv2.imread('noise_image.jpg', cv2.IMREAD_GRAYSCALE)

resized_image = cv2.resize(noisy_image, (512, 512))
# Normalize the image to [0, 1] range
resized_image_norm = resized_image.astype(np.float32) / 255.0

# Define function for filtering in image space (Gaussian blur)
def image_filtering(image):
    denoised_image = cv2.GaussianBlur(image, (5, 5), 0)
    return denoised_image

# Denoise using image filtering
denoised_image = image_filtering(resized_image_norm)

# Calculate MSE, PSNR, and SSIM
mse = mean_squared_error(resized_image_norm, denoised_image)
psnr = peak_signal_noise_ratio(resized_image_norm, denoised_image, data_range=resized_image_norm.max() - resized_image_norm.min())
ssim = structural_similarity(resized_image_norm, denoised_image, data_range=resized_image_norm.max() - resized_image_norm.min())

print(f"Image Filtering: MSE={mse}, PSNR={psnr}, SSIM={ssim}")
# Compute the 2D Fourier Transform
ft_image = np.fft.fft2(resized_image)
ft_image_shifted = np.fft.fftshift(ft_image)

# Display the magnitude spectrum of the Fourier transform
magnitude_spectrum = 20 * np.log(np.abs(ft_image_shifted))
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.imshow(resized_image, cmap='gray')
plt.title('Original Image')
plt.subplot(1, 2, 2)
plt.imshow(magnitude_spectrum, cmap='gray')
plt.title('Magnitude Spectrum')
plt.colorbar()
plt.show()

# Filtering in Frequency Space with different low-pass filters
cutoffs = [20, 40, 60, 80, 100]  # Different cutoff frequencies for low-pass filters
plt.figure(figsize=(15, 5))
plt.subplot(1, len(cutoffs) + 1, 1)
plt.imshow(resized_image, cmap='gray')
plt.title('Original Image')
for i, cutoff in enumerate(cutoffs):
    rows, cols = resized_image.shape
    crow, ccol = rows // 2, cols // 2

    # Create a low-pass filter
    mask = np.zeros((rows, cols), np.uint8)
    mask[crow - cutoff:crow + cutoff, ccol - cutoff:ccol + cutoff] = 1

    # Apply the filter to the shifted Fourier transform
    ft_image_filtered = ft_image_shifted * mask

    # Inverse Fourier Transform
    ift_image = np.fft.ifftshift(ft_image_filtered)
    denoised_image = np.fft.ifft2(ift_image)
    denoised_image = np.abs(denoised_image)

    # Display the denoised image
    plt.subplot(1, len(cutoffs) + 1, i + 2)
    plt.imshow(denoised_image, cmap='gray')
    plt.title(f'Cutoff: {cutoff}')
plt.show()
