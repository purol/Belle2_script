import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([1.50833333333, 2.62500000001, 3.74166666668])
x_error = np.array([0.55833333333, 0.55833333333, 0.55833333333])
y1 = np.array([0.000025, 2.319038, 0.00001])  # data
y1_error = np.array([0.018579, 2.075404, 0.010699])
y2 = np.array([3.019785, 4.636813, 1.685280])  # MC (corrected)
y2_error = np.array([0.736193, 0.902065, 0.554836])

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
