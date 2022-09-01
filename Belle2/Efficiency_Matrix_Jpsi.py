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
gen_list = [r"$B^{+}\rightarrow K^{+}J/\psi$",
            r"$B^{+}\rightarrow K^{+}\pi^{0}J/\psi$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}J/\psi$",      
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}J/\psi$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{0}J/\psi$",
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}J/\psi$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{+}J/\psi$",   
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{+}J/\psi$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{+}\pi^{0}J/\psi$",
            r"$B^{+}\rightarrow K^{+}\pi^{0}\pi^{0}J/\psi$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{0}\pi^{0}J/\psi$",
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}\pi^{0}J/\psi$",
            r"$B^{+}\rightarrow K^{+}K^{-}K^{+}J/\psi$",
            r"$B^{+}\rightarrow K^{+}K^{-}K^{0}\pi^{+}J/\psi$",
            r"$B^{+}\rightarrow K^{+}K^{-}K^{+}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}J/\psi$",
            r"$B^{0}\rightarrow K^{0}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}J/\psi$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}J/\psi$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{+}\pi^{-}J/\psi$",
            r"$B^{0}\rightarrow K^{0}\pi^{0}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{0}\pi^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{+}K^{-}K^{0}J/\psi$",
            r"$B^{0}\rightarrow K^{+}K^{-}K^{+}\pi^{-}J/\psi$",
            r"$B^{0}\rightarrow K^{+}K^{-}K^{0}\pi^{0}J/\psi$",
            r"$others$"
            ]

correction = []
correction.append(1353401.0)
correction.append(880864.0)
correction.append(1753251.0)
correction.append(1973940.0)
correction.append(2214036.0)
correction.append(1223776.0)
correction.append(461945.0)
correction.append(286674.0)
correction.append(378742.0)
correction.append(452343.0)
correction.append(321948.0)
correction.append(341842.0)
correction.append(77819.0)
correction.append(42265.0)
correction.append(21917.0)
correction.append(1.0)
correction.append(1659619.0)
correction.append(295577.0)
correction.append(1971422.0)
correction.append(2604225.0)
correction.append(960194.0)
correction.append(1155456.0)
correction.append(363364.0)
correction.append(170813.0)
correction.append(251427.0)
correction.append(306204.0)
correction.append(270138.0)
correction.append(61690.0)
correction.append(28477.0)
correction.append(21410.0)
correction.append(7591203.0)

for i in range(0,30):
    for j in range(0,31):
        values[i][j] = (values[i][j] / correction[i])

df_cm = pd.DataFrame(values, index = [i for i in reco_list], columns = [i for i in gen_list])
#df_cm = pd.DataFrame(confusion_matrix)

plt.figure(figsize = (36,24))
ax = sn.heatmap(df_cm, annot=True, cmap="YlGnBu", fmt='.1e')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize = 18, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), fontsize = 18)
plt.savefig("eff_matrix_Jpsi.png", bbox_inches='tight')
