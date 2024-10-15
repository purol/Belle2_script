#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt


confusion_matrix = pd.read_excel("./data.xlsx", header=None, index_col=None, nrows=3)
#print(confusion_matrix)
confusion_matrix.to_numpy()

#axis_list = [r"$\tau^{\mp}\rightarrow\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\pi^{0}\nu$",r"$Z\rightarrow q \bar{q}$"]
reco_list = [r"$0.0 < M_{X_{s}}^{reco} < 0.6 GeV$",
             r"$0.6 < M_{X_{s}}^{reco} < 1.0 GeV$",
             r"$1.0 < M_{X_{s}}^{reco} < 2.0 GeV$"
             ]
gen_list = [r"$0.0 < M_{X_{s}}^{true} < 0.6 GeV$",
            r"$0.6 < M_{X_{s}}^{true} < 1.0 GeV$",
            r"$1.0 GeV < M_{X_{s}}^{true}$"
            ]


df_cm = pd.DataFrame(confusion_matrix.values, index = [i for i in reco_list], columns = [i for i in gen_list])
#df_cm = pd.DataFrame(confusion_matrix)

plt.figure(figsize = (12,8))
ax = sn.heatmap(df_cm, annot=True, cmap="YlGnBu", annot_kws={'size': 25})
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize = 18, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), rotation=45, fontsize = 18)
plt.savefig("migration.png", bbox_inches='tight')
