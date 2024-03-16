import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([2.1125, 2.7375, 3.3625, 3.9875])
x_error = np.array([0.3125, 0.3125, 0.3125, 0.3125])
y1 = np.array([12.156990, 4.882687, 4.410669, 0])  # data
y1_error = np.array([5.980024, 4.319868, 3.763518, 0])
y2 = np.array([10.031623, 17.807606, 9.050026, 0.372727])  # MC (corrected)
y2_error = np.array([1.253486, 1.685581, 1.203908, 0.216617])

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
