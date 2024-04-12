import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
x_error = np.array([0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000])
y1 = np.array([0.000000, 0.000000, 0.000000, 0.005084, 0.002757, 0.000000, 0.000000, 0.000000, 0.000000, 0.003316, 0.000000, 0.000000, 0.000000, 0.004030, 0.000000, 0.000000, 0.010556, 0.000000, 0.000869, 0.002283, 0.000567, 0.002739, 0.009615, 0.004677])  # data
y1_error = np.array([0.000000, 0.000000, 0.000000, 0.003629, 0.002761, 0.000000, 0.000000, 0.000000, 0.000000, 0.003321, 0.000000, 0.000000, 0.000000, 0.004040, 0.000000, 0.000000, 0.010590, 0.000000, 0.000870, 0.001618, 0.000567, 0.001979, 0.004622, 0.004689])

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
plt.xlim(-0.5, 23.5)

plt.tight_layout()
plt.show()
