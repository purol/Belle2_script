#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

# KID
file_KID_cov_all = "KID_cov_all.txt"
file_KID_cov_partial = "KID_cov_partial.txt"

KID_cov_all = np.loadtxt(file_KID_cov_all, dtype=float)
KID_cov_partial = np.loadtxt(file_KID_cov_partial, dtype=float)

ar = np.array(KID_cov_all)
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-0.00002, vmax=0.0003)
ax.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('KID_cov_all.png')

plt.clf()

ar_p = np.array(KID_cov_partial)
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[], vmin=-0.00002, vmax=0.0003)
ax_p.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('KID_cov_partial.png')

plt.clf()

# PID
file_PID_cov_all = "PID_cov_all.txt"
file_PID_cov_partial = "PID_cov_partial.txt"

PID_cov_all = np.loadtxt(file_PID_cov_all, dtype=float)
PID_cov_partial = np.loadtxt(file_PID_cov_partial, dtype=float)

ar = np.array(PID_cov_all)
ax = sn.heatmap(ar, cmap="viridis", xticklabels=[], yticklabels=[])#,vmin=-0.00002, vmax=0.0003)
ax.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('PID_cov_all.png')

plt.clf()

ar_p = np.array(PID_cov_partial)
ax_p = sn.heatmap(ar_p, cmap="viridis", xticklabels=[], yticklabels=[])#,vmin=-0.00002, vmax=0.0003)
ax_p.tick_params(left=False, bottom=False)
#plt.show()
plt.savefig('PID_cov_partial.png')
