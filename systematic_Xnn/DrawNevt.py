import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([2.1125, 2.7375, 3.3625, 3.9875])
x_error = np.array([0.3125, 0.3125, 0.3125, 0.3125])
y1 = np.array([10.172406, 5.009899, 3.125226, 0])  # data
y1_error = np.array([5.854116, 3.948315, 3.334863, 0])
y2 = np.array([10.134165, 17.989634, 9.144640, 0.376537])  # MC (corrected)
y2_error = np.array([1.266299, 1.702811, 1.216369, 0.218832])

# Define an offset for x-values
offset = 0.0

# Plotting the histograms as points with offset and x-error bars
plt.errorbar(x - offset, y1, xerr=x_error, yerr=y1_error, fmt='bo', label='Data', capsize=5, elinewidth=2, alpha=0.7)
plt.errorbar(x + offset, y2, xerr=x_error, yerr=y2_error, fmt='ro', label='signal MC with correction', capsize=5, elinewidth=2, alpha=0.7)

plt.title(r'$B \rightarrow X_{s} p \bar{p}$ analysis')
plt.xlabel(r'$M_{p\bar{p}}$ [GeV]')
plt.ylabel('the number of events')
plt.legend()

# Adding a grid
plt.grid(True)

# Setting x-axis limits
plt.xlim(1.8, 4.3)

plt.tight_layout()
plt.show()
