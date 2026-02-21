import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
x_error = np.array([0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000])
y1 = np.array([0.002333, 0.001047, 0.003299, 0.000000, 0.000000, 0.000000, 0.000000, 0.002058, 0.004117, 0.000000, 0.004250, 0.003194, 0.000386, 0.005024, 0.008507])  # data
y1_error = np.array([0.001099, 0.000609, 0.003303, 0.000000, 0.000000, 0.000000, 0.000000, 0.002061, 0.004122, 0.000000, 0.001252, 0.001149, 0.000386, 0.002316, 0.005087])

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
