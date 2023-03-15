#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

# pi0
file_pi0_cov_all = "pi0_corr_all.txt"
file_pi0_cov_partial = "pi0_corr_partial.txt"

pi0_cov_all = np.loadtxt(file_pi0_cov_all, dtype=float)
pi0_cov_partial = np.loadtxt(file_pi0_cov_partial, dtype=float)

ar = np.array(pi0_cov_all)
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-0.3, vmax=1.0)
ax.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('pi0_corr_all.png')

plt.clf()

ar_p = np.array(pi0_cov_partial)
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-0.3, vmax=1.0)
ax_p.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('pi0_corr_partial.png')

plt.clf()
