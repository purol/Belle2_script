import matplotlib.pyplot as plt
import numpy as np

# Input data
x_boundary = [1.15, 1.5, 2.0, 2.4, 3.0]
y1 = [0.001650, 0.001382, 0.000303, 0.000032] # PYTHIA
y1_error = [0.000042, 0.000044, 0.000044, 0.000032]
y2 = [0.001714, 0.001413, 0.000154, np.nan] #K1
y2_error = [0.000076, 0.000215, 0.000154, np.nan]

#y1 = [0.000522, 0.000495, 0.000283, 0.000052] # PYTHIA
#y1_error = [0.000022, 0.000023, 0.000038, 0.000037]
#y2 = [0.000618, 0.000484, np.nan, np.nan] #K1
#y2_error = [0.00004, 0.000099, np.nan, np.nan]

#y1 = [0.000225, 0.000203, 0.000057, np.nan] # PYTHIA
#y1_error = [0.000016, 0.000017, 0.000019, np.nan]
#y2 = [0.000295, 0.000166, np.nan, np.nan] #K1
#y2_error = [0.000031, 0.000074, np.nan, np.nan]

#y1 = [0.001389, 0.001219, 0.000493, 0.000513] # PYTHIA
#y1_error = [0.000035, 0.000035, 0.000050, 0.000115]
#y2 = [0.001439, 0.001299, 0.000455, np.nan] #K1
#y2_error = [0.00006, 0.000162, 0.000204, np.nan]

y1_legend = r"$B \rightarrow X_{s}\nu\bar{\nu}$ MC with Pythia"
y2_legend = r"$B \rightarrow K_{1}(1270)\nu\bar{\nu}$ MC"

# Compute bin centers and half-widths
x_centers = [(x_boundary[i] + x_boundary[i+1]) / 2 for i in range(len(x_boundary)-1)]
xerr = [(x_boundary[i+1] - x_boundary[i]) / 2 for i in range(len(x_boundary)-1)]

# Plot with both vertical and horizontal error bars
plt.figure(figsize=(8, 5))
plt.errorbar(x_centers, y1, xerr=xerr, yerr=y1_error, fmt='o', label=y1_legend, capsize=5)
plt.errorbar(x_centers, y2, xerr=xerr, yerr=y2_error, fmt='s', label=y2_legend, capsize=5)

# Optional: draw bin boundaries
for x in x_boundary:
    plt.axvline(x=x, color='gray', linestyle='--', linewidth=0.5)

# Labels and legend
plt.xlabel(r'$M_{X_{s}}^{true}$ [GeV]')
plt.ylabel('efficiency')
plt.title('')
plt.legend()
plt.tight_layout()
plt.show()
