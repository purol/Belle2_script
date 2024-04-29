#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import FixedLocator

Nbins = 12

# Fragmentation
file_Fragmentation_corr_all = "Fragmentation_corr_all.txt"
file_Fragmentation_corr_partial = "Fragmentation_corr_partial.txt"

Fragmentation_corr_all = np.loadtxt(file_Fragmentation_corr_all, dtype=float)
Fragmentation_corr_partial = np.loadtxt(file_Fragmentation_corr_partial, dtype=float)

ar = np.array(Fragmentation_corr_all)
ar = np.delete(ar, np.s_[:Nbins * 2], axis = 0) # select SIGNAL only
ar = np.delete(ar, np.s_[:Nbins * 2], axis = 1) # select SIGNAL only
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax.tick_params(left=False, bottom=False)

labels = ["Signal"]
ticks = [ i * Nbins for i in range(2)]
ax.yaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax.set_yticklabels(labels, rotation=90, va="center")
ax.xaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax.set_xticklabels(labels, rotation=0)

#plt.show()
plt.savefig('Fragmentation_corr_all.png')

plt.clf()

ar_p = np.array(Fragmentation_corr_partial)
ar_p = np.delete(ar_p, np.s_[:Nbins * 2], axis = 0) # select SIGNAL only
ar_p = np.delete(ar_p, np.s_[:Nbins * 2], axis = 1) # select SIGNAL only
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-1.0, vmax=1.0)
ax_p.tick_params(left=False, bottom=False)

labels = ["Signal"]
ticks = [i * Nbins for i in range(2)]
ax_p.yaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax_p.set_yticklabels(labels, rotation=90, va="center")
ax_p.xaxis.set_major_locator(FixedLocator([(t0 + t1) / 2 for t0, t1 in zip(ticks[:-1], ticks[1:])]))
ax_p.set_xticklabels(labels, rotation=0)

#plt.show()
plt.savefig('Fragmentation_corr_partial.png')

plt.clf()
