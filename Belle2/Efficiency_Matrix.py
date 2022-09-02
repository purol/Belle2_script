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
gen_list = [r"$B^{+}\rightarrow K^{+}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{+}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\nu\bar{\nu}$",      
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",   
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{+}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{0}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{+}K^{-}K^{+}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{+}K^{-}K^{0}\pi^{+}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{+}K^{-}K^{+}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{0}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{+}\pi^{-}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{0}\pi^{0}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\pi^{0}\pi^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}K^{-}K^{0}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}K^{-}K^{+}\pi^{-}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{+}K^{-}K^{0}\pi^{0}\nu\bar{\nu}$",
            r"$others$"
            ]

correction = []
correction.append(1350.0)
correction.append(1651.6)
correction.append(3300.4)
correction.append(673.9)
correction.append(825.3)
correction.append(473.3)
correction.append(119.1)
correction.append(24.6)
correction.append(148.2)
correction.append(127.6)
correction.append(102.9)
correction.append(91.8)
correction.append(54.8)
correction.append(21.7)
correction.append(11.0)
correction.append(1184.0)
correction.append(3034.7)
correction.append(1515.9)
correction.append(835.5)
correction.append(677.4)
correction.append(120.5)
correction.append(475.9)
correction.append(149.8)
correction.append(25.0)
correction.append(129.2)
correction.append(104.4)
correction.append(92.4)
correction.append(54.9)
correction.append(14.5)
correction.append(10.9)
correction.append(1738.5)

for i in range(0,30):
    for j in range(0,31):
        values[i][j] = (values[i][j] / correction[i])

df_cm = pd.DataFrame(values, index = [i for i in reco_list], columns = [i for i in gen_list])
#df_cm = pd.DataFrame(confusion_matrix)

plt.figure(figsize = (36,24))
ax = sn.heatmap(df_cm, annot=True, cmap="YlGnBu", fmt='.1e')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize = 18, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), fontsize = 18)
plt.xlabel(r"MC decay modes", fontsize=25)
plt.ylabel(r"Reconstructed decay modes", fontsize=25)
plt.savefig("eff_matrix.png", bbox_inches='tight')
