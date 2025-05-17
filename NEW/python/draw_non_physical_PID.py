import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
x_error = np.array([0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000])
y1 = np.array([0.002371, 0.000924, 0.003084, 0.000000, 0.000000, 0.000000, 0.000000, 0.001973, 0.004319, 0.000000, 0.003971, 0.002977, 0.000360, 0.004579, 0.008279])  # data
y1_error = np.array([0.001123, 0.000537, 0.003089, 0.000000, 0.000000, 0.000000, 0.000000, 0.001976, 0.004325, 0.000000, 0.001169, 0.001070, 0.000360, 0.002116, 0.005053])

# Define an offset for x-values
offset = 0.0

# Plotting the histograms as points with offset and x-error bars
plt.errorbar(x, y1, xerr=x_error, yerr=y1_error, fmt='bo', label='', capsize=5, elinewidth=2, alpha=0.7)

plt.title(r'ratio of event with non-physical PID correction')
plt.xlabel(r'index')
plt.ylabel('the ratio of event')

# Adding a grid
plt.grid(True)

# Setting x-axis limits
plt.xlim(-0.5, 14.5)

plt.tight_layout()
plt.show()
