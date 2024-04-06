import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0])
x_error = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
y1 = np.array([0.002892, 0.0, 0.0, 0.003056, 0.0, 0.0, 0.0, 0.017146, 0.0, 0.0, 0.0, 0.0, 0.004952, 0.002605, 0.015318, 0.003256, 0.0, 0.0])  # data
y1_error = np.array([0.002897, 0.0, 0.0, 0.003061, 0.0, 0.0, 0.0, 0.017235, 0.0, 0.0, 0.0, 0.0, 0.003553, 0.002610, 0.008133, 0.003262, 0.0, 0.0])

# Define an offset for x-values
offset = 0.0

# Plotting the histograms as points with offset and x-error bars
plt.errorbar(x, y1, xerr=x_error, yerr=y1_error, fmt='bo', label='', capsize=5, elinewidth=2, alpha=0.7)

plt.title(r'ratio of event with non-physical PID correction')
plt.xlabel(r'bin number')
plt.ylabel('the ratio of event')
plt.legend()

# Adding a grid
plt.grid(True)

# Setting x-axis limits
plt.xlim(-0.5, 17.5)

plt.tight_layout()
plt.show()
