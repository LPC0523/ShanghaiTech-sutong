import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tifffile

# Load the TIFF image
image = tifffile.imread('MTF-slice.tif')

# Identify the center of the wire
center_x = image.shape[1] // 2
center_y = image.shape[0] // 2

# Draw a line passing through the center of the wire
line_points = np.linspace(0, image.shape[1]-1, image.shape[1]).astype(int)

# Extract intensity values along this line
intensity_profile = [image[center_y, x] for x in line_points]

# Subtract background
background = np.min(intensity_profile)
intensity_profile -= background

# Plot the intensity profile
plt.plot(line_points, intensity_profile, label='Intensity Profile')
plt.xlabel('Position along Line')
plt.ylabel('Intensity')
plt.title('Line Intensity Profile of the Wire')
plt.grid(True)
plt.show()

# Fit the intensity profile with a Gaussian function
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Initial guess for parameters
initial_guess = [np.max(intensity_profile), np.argmax(intensity_profile), 10]

# Perform curve fitting
popt, pcov = curve_fit(gaussian, line_points, intensity_profile, p0=initial_guess)

# Plot the fitted Gaussian
plt.plot(line_points, intensity_profile, label='Intensity Profile')
plt.plot(line_points, gaussian(line_points, *popt), 'r--', label='Gaussian Fit')
plt.xlabel('Position along Line')
plt.ylabel('Intensity')
plt.title('Line Intensity Profile of the Wire with Gaussian Fit')
plt.grid(True)
plt.legend()
plt.show()

# Print fitted parameters
print("Fitted parameters:")
print("Amplitude (A):", popt[0])
print("Mean (mu):", popt[1])
print("Standard deviation (sigma):", popt[2])

# Calculate Line Spread Function (LSF)
lsf = gaussian(line_points, *popt)

# Calculate Modulation Transfer Function (MTF) by taking the Fourier Transform of LSF
mtf = np.fft.fft(lsf)

# Normalize MTF
mtf /= np.abs(mtf).max()

# Spatial frequencies
frequencies = np.fft.fftfreq(len(mtf), d=line_points[1] - line_points[0])

# Plot the MTF
plt.plot(frequencies, np.abs(mtf), label='MTF')
plt.xlabel('Spatial Frequency')
plt.ylabel('Modulation Transfer Function (MTF)')
plt.title('Modulation Transfer Function (MTF)')
plt.grid(True)
plt.show()

# Fit the MTF with a Gaussian function
def gaussian_mtf(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Initial guess for parameters
initial_guess_mtf = [np.max(np.abs(mtf)), 0, 0.1]  # Use peak value of MTF as initial guess for amplitude

# Perform curve fitting
popt_mtf, pcov_mtf = curve_fit(gaussian_mtf, frequencies, np.abs(mtf), p0=initial_guess_mtf)

# Plot the fitted Gaussian
plt.plot(frequencies, np.abs(mtf), label='MTF')
plt.plot(frequencies, gaussian_mtf(frequencies, *popt_mtf), 'r--', label='Gaussian Fit')
plt.xlabel('Spatial Frequency')
plt.ylabel('Modulation Transfer Function (MTF)')
plt.title('Modulation Transfer Function (MTF) with Gaussian Fit')
plt.grid(True)
plt.legend()
plt.show()

# Calculate spatial resolution at 10% MTF
ten_percent_mtf = 0.1 * popt_mtf[0]  # 10% of the peak value
ten_percent_mtf_frequency = frequencies[np.argmin(np.abs(np.abs(mtf) - ten_percent_mtf))]
spatial_resolution = 1 / (2 * ten_percent_mtf_frequency)

# Print spatial resolution at 10% MTF
print("Spatial resolution at 10% MTF:", spatial_resolution)
