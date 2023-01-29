#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

# BR
file_BR_corr_all = "BR_corr_all.txt"
file_BR_corr_partial = "BR_corr_partial.txt"

BR_corr_all = np.loadtxt(file_BR_corr_all, dtype=float)
BR_corr_partial = np.loadtxt(file_BR_corr_partial, dtype=float)

ar = np.array(BR_corr_all)
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[])#, vmin=-0.3, vmax=1.0)
ax.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('BR_corr_all.png')

plt.clf()

ar_p = np.array(BR_corr_partial)
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[])#, vmin=0.0, vmax=0.0001)
ax_p.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('BR_corr_partial.png')

plt.clf()
