import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
x_error = np.array([0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000])
y1 = np.array([0.001852, 0.001361, 0.003112, 0.000000, 0.000000, 0.000000, 0.000000, 0.002003, 0.004400, 0.000000, 0.003059, 0.002560, 0.001243, 0.004635, 0.008051])  # data
y1_error = np.array([0.000995, 0.000694, 0.003116, 0.000000, 0.000000, 0.000000, 0.000000, 0.002005, 0.004406, 0.000000, 0.001052, 0.000979, 0.000953, 0.002140, 0.004927])

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
