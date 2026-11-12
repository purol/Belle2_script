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
correction.append(1995.975877)
correction.append(1502.585535)
correction.append(2988.635366)
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
correction.append(1748.551733)
correction.append(2698.236699)
correction.append(1344.974956)
correction.append(1395.602934)
correction.append(861.034692)
correction.append(280.510287)
correction.append(751.723680)
correction.append(131.773694)
correction.append(25.787289)
correction.append(183.495400)
correction.append(179.368498)
correction.append(142.182852)
correction.append(44.223493)
correction.append(13.812601)
correction.append(10.466336)
correction.append(1943.409873)

Efficiency_diagonal = []
Efficiency = np.zeros((30,30))
cross_feed = np.zeros((30,30))

for i in range(0,30):
    Efficiency_diagonal.append(values[i][i] / correction[i])

for i in range(0,30):
    for j in range(0,30):
        Efficiency[i][j] = (values[i][j] / correction[j])

for i in range(0,30):
    for j in range(0,30):
        if Efficiency_diagonal[j] != 0:
            cross_feed[i][j] = (Efficiency[i][j] / Efficiency_diagonal[j])
        else:
            cross_feed[i][j] = np.nan


print(Efficiency_diagonal[0])


reco_list.reverse()
Efficiency_diagonal.reverse()

plt.rcParams['figure.figsize'] = [15, 12]

fig, ax = plt.subplots()
bars = ax.barh(reco_list, Efficiency_diagonal)
#ax.bar_label(bars)

#plt.barh(y, Nevt)
#plt.yticks(y, decay)

plt.xlabel("Efficiency", fontsize=20)
plt.ylabel("Decay Modes", fontsize=20)
plt.grid(True, axis = "x", linestyle="--")

ax.tick_params(axis='x', labelsize=16)  # Change font size of x-axis tick labels
ax.tick_params(axis='y', labelsize=16)  # Change font size of y-axis tick labels

plt.savefig("Diagonal_Efficiency.pdf", bbox_inches='tight')

# draw cross-feed matrix
reco_list.reverse()

df_cm = pd.DataFrame(cross_feed, index = [i for i in reco_list], columns = [i for i in reco_list])

with plt.rc_context({
    "pdf.fonttype": 42,
    "mathtext.fontset": "stix",
}):
    # Match the final width: 0.98 * 15.92 cm in the thesis.
    fig_width_cm = 15.6016
    fig_height_cm = 15.0
    # Equal physical width and height for the 30-by-30 matrix.
    matrix_height = 0.72 * fig_width_cm / fig_height_cm
    fig = plt.figure(
        figsize=(fig_width_cm / 2.54, fig_height_cm / 2.54),
        dpi=300,
    )

    ax = fig.add_axes([0.21, 0.23, 0.72, matrix_height])
    cax = fig.add_axes([0.945, 0.23, 0.018, matrix_height])

    sn.heatmap(
        df_cm,
        ax=ax,
        cbar_ax=cax,
        annot=True,
        # Preserve the approximately 2.2 pt cell numbers on the thesis page.
        annot_kws={"fontsize": 2.2},
        cmap="YlGnBu",
        fmt=".1e",
        square=True,
        xticklabels=True,
        yticklabels=True,
    )
    ax.set_facecolor("lightgray")

    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=60,
        rotation_mode="anchor",
        fontsize=7.5,
        ha="right",
        va="top",
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=7.5,
    )

    ax.set_xlabel("MC decay modes", fontsize=9, labelpad=4)
    ax.set_ylabel(
        "Reconstructed decay modes",
        fontsize=9,
        labelpad=4,
    )

    ax.tick_params(axis="both", length=1.5, width=0.5, pad=2)
    cax.tick_params(labelsize=8.2, length=2, width=0.5, pad=2)

    fig.savefig("cross_feed_matrix.pdf", bbox_inches=None)
    plt.close(fig)
