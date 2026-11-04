#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import FixedLocator

# Fonts are enlarged to remain readable when scaled to thesis column width.
plt.rcParams.update({
    "font.size": 26,
    "xtick.labelsize": 26,
    "ytick.labelsize": 26,
    "axes.labelsize": 26,
    "axes.titlesize": 26,
})

plt.figure(figsize=(6.4 * 1.6, 4.8 * 1.6))
Nbins = 15

# Fragmentation
file_Fragmentation_corr_all = "Fragmentation_corr_all.txt"
file_Fragmentation_corr_partial = "Fragmentation_corr_partial.txt"

Fragmentation_corr_all = np.loadtxt(file_Fragmentation_corr_all, dtype=float)
Fragmentation_corr_partial = np.loadtxt(file_Fragmentation_corr_partial, dtype=float)

ar = np.array(Fragmentation_corr_all)
ar = np.delete(ar, np.s_[:Nbins * 4], axis = 0) # select SIGNAL only
ar = np.delete(ar, np.s_[:Nbins * 4], axis = 1) # select SIGNAL only
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax.tick_params(left=False, bottom=False)

# Show selected bin indices at the centers of their cells.
bin_labels = np.array([1, 5, 10, 15])
positions = bin_labels - 0.5
ax.set_xticks(positions)
ax.set_yticks(positions)
ax.set_xticklabels(bin_labels, rotation=0)
ax.set_yticklabels(bin_labels, rotation=0)
ax.set_xlabel("Bin index")
ax.set_ylabel("Bin index")
ax.set_title(r"$S_{3}$")

#plt.show()
cbar = ax.collections[0].colorbar
cbar.set_label("Correlation coefficient")
cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
plt.tight_layout()
plt.savefig('Fragmentation_corr_all.pdf', dpi=300, bbox_inches="tight")

plt.clf()

ar_p = np.array(Fragmentation_corr_partial)
ar_p = np.delete(ar_p, np.s_[:Nbins * 4], axis = 0) # select SIGNAL only
ar_p = np.delete(ar_p, np.s_[:Nbins * 4], axis = 1) # select SIGNAL only
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax_p.tick_params(left=False, bottom=False)

# Show selected bin indices at the centers of their cells.
bin_labels = np.array([1, 5, 10, 15])
positions = bin_labels - 0.5
ax_p.set_xticks(positions)
ax_p.set_yticks(positions)
ax_p.set_xticklabels(bin_labels, rotation=0)
ax_p.set_yticklabels(bin_labels, rotation=0)
ax_p.set_xlabel("Bin index")
ax_p.set_ylabel("Bin index")
ax_p.set_title(r"$S_{3}$")

#plt.show()
cbar = ax_p.collections[0].colorbar
cbar.set_label("Correlation coefficient")
cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
plt.tight_layout()
plt.savefig('Fragmentation_corr_partial.pdf', dpi=300, bbox_inches="tight")

plt.clf()
