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
correction.append(1368900.0)
correction.append(890767.0)
correction.append(1773288.0)
correction.append(1996210.0)
correction.append(2238907.0)
correction.append(1237632.0)
correction.append(467178.0)
correction.append(289966.0)
correction.append(383141.0)
correction.append(457464.0)
correction.append(325602.0)
correction.append(345567.0)
correction.append(78717.0)
correction.append(42725.0)
correction.append(22168.0)
correction.append(1.0)
correction.append(1678648.0)
correction.append(298908.0)
correction.append(1993814.0)
correction.append(2633753.0)
correction.append(971250.0)
correction.append(1168624.0)
correction.append(367522.0)
correction.append(172739.0)
correction.append(254144.0)
correction.append(309677.0)
correction.append(273233.0)
correction.append(62368.0)
correction.append(28790.0)
correction.append(21625.0)
correction.append(7677320.0)

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
