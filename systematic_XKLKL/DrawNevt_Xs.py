import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([2.625])
x_error = np.array([1.675])
y1 = np.array([1.414669])  # data
y1_error = np.array([3.746412])
y2 = np.array([9.854887])  # MC (corrected)
y2_error = np.array([1.340821])

# Define an offset for x-values
offset = 0.0

# Plotting the histograms as points with offset and x-error bars
plt.errorbar(x - offset, y1, xerr=x_error, yerr=y1_error, fmt='bo', label='Data', capsize=5, elinewidth=2, alpha=0.7)
plt.errorbar(x + offset, y2, xerr=x_error, yerr=y2_error, fmt='ro', label='signal MC with correction', capsize=5, elinewidth=2, alpha=0.7)

plt.title(r'$B \rightarrow X_{s} K_{S}^{0} K_{S}^{0}$ analysis')
plt.xlabel(r'$M_{K_{S}^{0}K_{S}^{0}}$ [GeV]')
plt.ylabel('the number of events')
plt.legend()

# Adding a grid
plt.grid(True)

# Setting x-axis limits
plt.xlim(0.95, 4.3)

plt.tight_layout()
plt.show()
