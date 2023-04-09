#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt


confusion_matrix = pd.read_excel("./efficiency_data.xlsx", header=None, index_col=None, nrows=30)
#print(confusion_matrix)
values = confusion_matrix.to_numpy(dtype='float', copy=True)

#axis_list = [r"$\tau^{\mp}\rightarrow\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\pi^{0}\nu$",r"$Z\rightarrow q \bar{q}$"]
reco_list = [r"$B^{+}\rightarrow K^{+}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}\pi^{0}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{0}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{+}\pi^{0}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}K^{-}K^{+}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}K^{-}K^{0}_{S}\pi^{+}\nu\bar{\nu}$",
             r"$B^{+}\rightarrow K^{+}K^{-}K^{+}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{0}_{S}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}\pi^{-}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{0}_{S}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{+}\pi^{-}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{0}_{S}\pi^{0}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{0}\pi^{0}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}K^{-}K^{0}_{S}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}K^{-}K^{+}\pi^{-}\nu\bar{\nu}$",
             r"$B^{0}\rightarrow K^{+}K^{-}K^{0}_{S}\pi^{0}\nu\bar{\nu}$"
             ]

correction = []
correction.append(2007.3)
correction.append(1893.5)
correction.append(3783.7)
correction.append(734.9)
correction.append(899.9)
correction.append(516.1)
correction.append(129.8)
correction.append(26.8)
correction.append(161.6)
correction.append(139.1)
correction.append(112.3)
correction.append(100.1)
correction.append(59.8)
correction.append(23.6)
correction.append(12.0)
correction.append(1760.2)
correction.append(3482.8)
correction.append(1739.7)
correction.append(920.1)
correction.append(746.0)
correction.append(132.7)
correction.append(524.0)
correction.append(164.9)
correction.append(27.5)
correction.append(142.3)
correction.append(114.9)
correction.append(101.8)
correction.append(60.4)
correction.append(16.0)
correction.append(12.0)
correction.append(1905.9)

Efficiency_diagonal = []
Efficiency = np.zeros((30,30))
cross_feed = np.zeros((30,30))

for i in range(0,30):
    Efficiency_diagonal.append(values[i][i] / correction[i])

for i in range(0,30):
    for j in range(0,30):
        Efficiency[i][j] = (values[i][j] / correction[i])

for i in range(0,30):
    for j in range(0,30):
        if Efficiency_diagonal[j] != 0:
            cross_feed[i][j] = (Efficiency[i][j] / Efficiency_diagonal[j])
        else:
            cross_feed[i][j] = 0

print(values[22][22])
print(values[20][20])
print(values[20][22])
print(Efficiency[20][22])
print(Efficiency_diagonal[22])
print(cross_feed[20][22])
print(Efficiency_diagonal[0])

reco_list.reverse()
Efficiency_diagonal.reverse()

plt.rcParams['figure.figsize'] = [15, 12]

fig, ax = plt.subplots()
bars = ax.barh(reco_list, Efficiency_diagonal)
#ax.bar_label(bars)

#plt.barh(y, Nevt)
#plt.yticks(y, decay)

plt.xlabel("Efficiency")
plt.ylabel("Decay Modes")
plt.grid(True, axis = "x", linestyle="--")

plt.savefig("Diagonal_Efficiency.png", bbox_inches='tight')

# draw cross-feed matrix
reco_list.reverse()

df_cm = pd.DataFrame(cross_feed, index = [i for i in reco_list], columns = [i for i in reco_list])

plt.figure(figsize = (36,24))
ax = sn.heatmap(df_cm, annot=True, cmap="YlGnBu", fmt='.1e')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize = 18, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), fontsize = 18)
plt.xlabel(r"MC decay modes", fontsize=25)
plt.ylabel(r"Reconstructed decay modes", fontsize=25)
plt.savefig("cross_feed_matrix.png", bbox_inches='tight')
