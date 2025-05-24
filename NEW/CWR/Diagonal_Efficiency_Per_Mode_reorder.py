#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

# reorder...
manual_order = [
    r"$B^{+}\rightarrow K^{+}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{0}_{S}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}\pi^{-}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{0}_{S}\pi^{0}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{0}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{0}_{S}\pi^{0}\pi^{0}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\pi^{0}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{+}\pi^{0}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{+}\pi^{-}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{0}_{S}\pi^{+}\pi^{-}\pi^{0}\pi^{0}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}K^{-}K^{+}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}K^{-}K^{0}_{S}\pi^{+}\nu\bar{\nu}$",
    r"$B^{+}\rightarrow K^{+}K^{-}K^{+}\pi^{0}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}K^{-}K^{0}_{S}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}K^{-}K^{+}\pi^{-}\nu\bar{\nu}$",
    r"$B^{0}\rightarrow K^{+}K^{-}K^{0}_{S}\pi^{0}\nu\bar{\nu}$"
]

confusion_matrix = pd.read_excel("./efficiency_data_MXs_larger_1.1.xlsx", header=None, index_col=None, nrows=30)
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
correction.append(0.002057)
correction.append(272.381492)
correction.append(539.575093)
correction.append(844.979127)
correction.append(1346.461865)
correction.append(737.196171)
correction.append(272.574339)
correction.append(25.318259)
correction.append(128.540240)
correction.append(180.494993)
correction.append(174.151800)
correction.append(140.068766)
correction.append(43.920483)
correction.append(20.925846)
correction.append(10.417798)
correction.append(0.001733)
correction.append(520.105275)
correction.append(260.003009)
correction.append(1395.602934)
correction.append(861.034692)
correction.append(280.509148)
correction.append(751.722541)
correction.append(131.773694)
correction.append(25.787289)
correction.append(183.495400)
correction.append(179.368498)
correction.append(142.182852)
correction.append(44.223493)
correction.append(13.811462)
correction.append(10.466336)
correction.append(1928.805056)

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


print(Efficiency_diagonal[0])

reco_list.reverse()
Efficiency_diagonal.reverse()

# reorder
index_map = {name: i for i, name in enumerate(reco_list)}
ordered_efficiencies = [Efficiency_diagonal[index_map[name]] for name in manual_order]

manual_order = manual_order[2:6] + ["  "] + manual_order[6:12] + ["   "] + manual_order[12:18] + ["    "] + manual_order[18:24] + ["     "] + manual_order[24:]
ordered_efficiencies = ordered_efficiencies[2:6] + [0] + ordered_efficiencies[6:12] + [0] + ordered_efficiencies[12:18] + [0] + ordered_efficiencies[18:24] + [0] + ordered_efficiencies[24:]

manual_order.reverse()
ordered_efficiencies.reverse()

plt.rcParams['figure.figsize'] = [15, 12]

fig, ax = plt.subplots()
bars = ax.barh(manual_order, ordered_efficiencies)
#ax.bar_label(bars)

#plt.barh(y, Nevt)
#plt.yticks(y, decay)

plt.xlabel("Efficiency", fontsize=20)
plt.ylabel("Decay Modes", fontsize=20)
plt.grid(True, axis = "x", linestyle="--")

ax.tick_params(axis='x', labelsize=16)  # Change font size of x-axis tick labels
ax.tick_params(axis='y', labelsize=16)  # Change font size of y-axis tick labels

ax.text(0.8, 0.05, r'$1.1$ GeV < $M_{X_{s}}$', transform=ax.transAxes, color='black', fontsize=16)

plt.savefig("Diagonal_Efficiency.png", bbox_inches='tight')
