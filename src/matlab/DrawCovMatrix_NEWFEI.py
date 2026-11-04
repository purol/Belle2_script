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

# pi0
file_pi0_cov_all = "NEWFEICAL_corr_all.txt"
file_pi0_cov_partial = "NEWFEICAL_corr_partial.txt"

pi0_cov_all = np.loadtxt(file_pi0_cov_all, dtype=float)
pi0_cov_partial = np.loadtxt(file_pi0_cov_partial, dtype=float)

ar = np.array(pi0_cov_all)
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax.tick_params(left=False, bottom=False)

labels = [
    r"$S_{1}$",
    r"$S_{2}$",
    r"$S_{3}$",
]
ticks = [i * Nbins for i in range(4)]
ax.yaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax.set_yticklabels(labels, rotation=0, va="center")
ax.xaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")

#plt.show()
cbar = ax.collections[0].colorbar
cbar.set_label("Correlation coefficient")
cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
plt.tight_layout()
plt.savefig('NEWFEICAL_corr_all.pdf', dpi=300, bbox_inches="tight")

plt.clf()

ar_p = np.array(pi0_cov_partial)
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax_p.tick_params(left=False, bottom=False)

labels = [
    r"$S_{1}$",
    r"$S_{2}$",
    r"$S_{3}$",
]
ticks = [i * Nbins for i in range(4)]
ax_p.yaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax_p.set_yticklabels(labels, rotation=0, va="center")
ax_p.xaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax_p.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")

#plt.show()
cbar = ax_p.collections[0].colorbar
cbar.set_label("Correlation coefficient")
cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
plt.tight_layout()
plt.savefig('NEWFEICAL_corr_partial.pdf', dpi=300, bbox_inches="tight")

plt.clf()
