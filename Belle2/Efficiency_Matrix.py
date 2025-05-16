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
correction.append(2001.326190
)
correction.append(1506.613291
)
correction.append(2996.646553
)
correction.append(847.244136
)
correction.append(1350.071124
)
correction.append(739.172262
)
correction.append(273.304988
)
correction.append(25.386126
)
correction.append(128.884799
)
correction.append(180.978819
)
correction.append(174.618623
)
correction.append(140.444228
)
correction.append(44.038214
)
correction.append(20.981938
)
correction.append(10.445723
)
correction.append(1752.974136
)
correction.append(2705.061027
)
correction.append(1348.376640
)
correction.append(1399.132666
)
correction.append(863.212404
)
correction.append(281.219748
)
correction.append(753.624925
)
correction.append(132.106974
)
correction.append(25.852509
)
correction.append(183.959493
)
correction.append(179.822153
)
correction.append(142.542458
)
correction.append(44.335342
)
correction.append(13.847536
)
correction.append(10.492807
)
correction.append(1948.325107)

for i in range(0,30):
    for j in range(0,31):
        values[i][j] = (values[i][j] / correction[i])

print(values[22][22])
print(values[20][20])
print(values[20][22])

df_cm = pd.DataFrame(values, index = [i for i in reco_list], columns = [i for i in gen_list])
#df_cm = pd.DataFrame(confusion_matrix)

plt.figure(figsize = (36,24))
ax = sn.heatmap(df_cm, annot=True, cmap="YlGnBu", fmt='.1e')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize = 18, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), fontsize = 18)
plt.xlabel(r"MC decay modes", fontsize=25)
plt.ylabel(r"Reconstructed decay modes", fontsize=25)
plt.savefig("eff_matrix.png", bbox_inches='tight')
