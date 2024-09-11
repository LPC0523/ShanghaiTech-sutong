import numpy as np
import matplotlib.pyplot as plt

# Define square wave function
def square_wave(x):
    if 0 <= x < 0.5:
        return 1
    elif 0.5 <= x < 1:
        return 0

# Define Fourier series coefficients for square wave
def a0():
    return 1

def an(n):
    if n % 2 != 0:
        return 2 / (np.pi * n)
    else:
        return 0

def bn(n):
    return 0

# Define function to reconstruct square wave using Fourier series
def reconstruct_square_wave(x, num_terms):
    result = a0() / 2
    for n in range(1, num_terms + 1):
        result += an(n) * np.cos(2 * np.pi * n * x) + bn(n) * np.sin(2 * np.pi * n * x)
    return result

# Generate x values
x_values = np.linspace(0, 2, 1000)

# Reconstruct square wave using Fourier series with different number of terms
num_terms_list = [1, 3, 5, 10]
plt.figure(figsize=(10, 6))
for num_terms in num_terms_list:
    reconstructed_wave = [reconstruct_square_wave(x, num_terms) for x in x_values]
    plt.plot(x_values, reconstructed_wave, label=f'{num_terms} Terms')

# Plot original square wave
original_wave = [square_wave(x) for x in x_values]
plt.plot(x_values, original_wave, label='Original Square Wave', linestyle='--', color='black')

plt.title('Reconstruction of Square Wave using Fourier Series')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
plt.show()
