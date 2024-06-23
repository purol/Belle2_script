#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import FixedLocator

plt.figure(figsize=(6.4 * 1.6, 4.8 * 1.6))
Nbins = 15

# BR
file_BR_corr_all = "FEI_corr_all.txt"
file_BR_corr_partial = "FEI_corr_partial.txt"

BR_corr_all = np.loadtxt(file_BR_corr_all, dtype=float)
BR_corr_partial = np.loadtxt(file_BR_corr_partial, dtype=float)

ar = np.array(BR_corr_all)
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax.tick_params(left=False, bottom=False)

labels = ["CHG", "MIX", "Signal_1", "Signal_2", "Signal_3"]
ticks = [i * Nbins for i in range(6)]
ax.yaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax.set_yticklabels(labels, rotation=90, va="center")
ax.xaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax.set_xticklabels(labels, rotation=0)

#plt.show()
plt.savefig('FEI_corr_all.png')

plt.clf()

ar_p = np.array(BR_corr_partial)
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax_p.tick_params(left=False, bottom=False)

labels = ["CHG", "MIX", "Signal_1", "Signal_2", "Signal_3"]
ticks = [i * Nbins for i in range(6)]
ax_p.yaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax_p.set_yticklabels(labels, rotation=90, va="center")
ax_p.xaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax_p.set_xticklabels(labels, rotation=0)

#plt.show()
plt.savefig('FEI_corr_partial.png')

plt.clf()
