import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.fftpack import fft, fftshift
from skimage import io

# Load the image
image = io.imread('MTF-slice.tif')

# Calculate the Line Intensity Profile (LIP)
center_y = image.shape[0] // 2
lip = image[center_y, :]

# Define Gaussian function for fitting
def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) / sigma)**2 / 2)

# Define range of x values
x = np.arange(len(lip))

# Fit the LIP with Gaussian function
popt, pcov = curve_fit(gaussian, x, lip, p0=[max(lip), len(lip) // 2, 1])

# Calculate FWHM (Full Width at Half Maximum) from sigma
fwhm = 2 * np.sqrt(2 * np.log(2)) * popt[2]

# Plot LIP and Gaussian fit
plt.figure(figsize=(10, 5))
plt.plot(x, lip, label='Line Intensity Profile (LIP)')
plt.plot(x, gaussian(x, *popt), 'r--', label='Gaussian Fit')
plt.xlabel('Position')
plt.ylabel('Intensity')
plt.title('Line Intensity Profile (LIP) and Gaussian Fit')
plt.legend()
plt.grid(True)
plt.show()

# Calculate Line Spread Function (LSF)
lsf = gaussian(x, *popt)

# Calculate Modulation Transfer Function (MTF)
mtf = np.abs(fftshift(fft(lsf)))
mtf /= mtf.max()  # Normalize MTF

# Find spatial resolution at 10% MTF
ten_percent_mtf = 0.1 * mtf.max()
spatial_resolution = np.argmax(mtf >= ten_percent_mtf)

# Plot MTF
plt.figure(figsize=(10, 5))
plt.plot(mtf, label='MTF')
plt.axhline(y=ten_percent_mtf, color='r', linestyle='--', label='10% MTF')
plt.xlabel('Spatial Frequency')
plt.ylabel('MTF')
plt.title('Modulation Transfer Function (MTF)')
plt.legend()
plt.grid(True)
plt.show()

print(f"Spatial resolution at 10% MTF: {spatial_resolution} pixels")
print(f"FWHM: {fwhm} pixels")
