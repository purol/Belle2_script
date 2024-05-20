import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
x_error = np.array([0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000])
y1 = np.array([0.002853, 0.000917, 0.003167, 0.000000, 0.000000, 0.000000, 0.000000, 0.002118, 0.004141, 0.000000, 0.003013, 0.002637, 0.001239, 0.003722, 0.009810])  # data
y1_error = np.array([0.001223, 0.000532, 0.003172, 0.000000, 0.000000, 0.000000, 0.000000, 0.002121, 0.004146, 0.000000, 0.001031, 0.001015, 0.000949, 0.001929, 0.005214])

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
