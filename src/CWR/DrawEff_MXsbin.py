import matplotlib.pyplot as plt
import numpy as np

# Input data
x_boundary = [1.15, 1.5, 2.0, 2.4, 3.0]
y1 = [0.001650, 0.001382, 0.000303, 0.000032] # PYTHIA
y1_error = [0.000042, 0.000044, 0.000044, 0.000032]
y2 = [0.001726, 0.001387, 0.000214, 0.000106] #K1
y2_error = [0.000024, 0.000067, 0.000057, 0.000075]
desc = r"$B^{+}\rightarrow K^{+} \pi^{-} \pi^{+} \nu \bar{\nu}$"

#y1 = [0.000522, 0.000495, 0.000283, 0.000052] # PYTHIA
#y1_error = [0.000022, 0.000023, 0.000038, 0.000037]
#y2 = [0.000598, 0.000522, 0.000192, 0.000219] #K1
#y2_error = [0.000012, 0.000033, 0.000042, 0.000083]
#desc = r"$B^{+}\rightarrow K^{0} \pi^{+} \pi^{0} \nu \bar{\nu}$"

#y1 = [0.000225, 0.000203, 0.000057, np.nan] # PYTHIA
#y1_error = [0.000016, 0.000017, 0.000019, np.nan]
#y2 = [0.000278, 0.000210, 0.000061, 0.000052] #K1
#y2_error = [0.000010, 0.000026, 0.000031, 0.000052]
#desc = r"$B^{0}\rightarrow K^{0} \pi^{+} \pi^{-} \nu \bar{\nu}$"

#y1 = [0.001389, 0.001219, 0.000493, 0.000513] # PYTHIA
#y1_error = [0.000035, 0.000035, 0.000050, 0.000115]
#y2 = [0.001474, 0.001225, 0.000498, 0.000186] #K1
#y2_error = [0.000019, 0.000050, 0.000068, 0.000076]
#desc = r"$B^{0}\rightarrow K^{+} \pi^{-} \pi^{0} \nu \bar{\nu}$"

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

plt.text(0.70, 0.80, desc, 
         transform=plt.gca().transAxes,  # use axes coords
         fontsize=12, 
         verticalalignment='top')

# Labels and legend
plt.xlabel(r'$M_{X_{s}}^{true}$ [GeV]')
plt.ylabel('efficiency')
plt.title('')
plt.legend()
plt.tight_layout()
plt.show()
