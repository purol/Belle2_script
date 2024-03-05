import matplotlib.pyplot as plt
import numpy as np

# Data
x = np.array([2.13333333333, 2.8, 3.46666666667])
x_error = np.array([1/3, 1/3, 1/3])
y1 = np.array([9.084567, 3.573505, 2.026306])  # data
y1_error = np.array([4.273256, 3.030385, 2.485834])
y2 = np.array([8.133702, 11.762898, 3.070198])  # MC (corrected)
y2_error = np.array([1.209195, 1.471909, 0.755383])
y3 = np.array([1.938013, 3.446400, 0.982748])  # peaking MC (corrected). CHG/MIX sample with Mbc > 5.27
y3_error = np.array([0.597488, 0.809191, 0.408252])

# Define an offset for x-values
offset = 0.0

# Plotting the histograms as points with offset and x-error bars
plt.errorbar(x - offset, y1, xerr=x_error, yerr=y1_error, fmt='bo', label='Data', capsize=5, elinewidth=2, alpha=0.7)
plt.errorbar(x         , y2, xerr=x_error, yerr=y2_error, fmt='ro', label='signal MC with correction', capsize=5, elinewidth=2, alpha=0.7)
plt.errorbar(x + offset, y3, xerr=x_error, yerr=y3_error, fmt='go', label='peaking MC with correction', capsize=5, elinewidth=2, alpha=0.7)

plt.title(r'$B \rightarrow X_{s} p \bar{p}$ analysis')
plt.xlabel(r'$M_{p\bar{p}}$ [GeV]')
plt.ylabel('the number of events')
plt.legend()

# Adding a grid
plt.grid(True)

# Setting x-axis limits
plt.xlim(1.8, 3.8)

plt.tight_layout()
plt.show()
