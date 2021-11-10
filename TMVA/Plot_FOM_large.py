#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt


confusion_matrix = pd.read_excel("./FOM_large.xlsx", header=None, index_col=None, nrows=20)
#print(confusion_matrix)
confusion_matrix.to_numpy()

#axis_list = [r"$\tau^{\mp}\rightarrow\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\pi^{0}\nu$",r"$Z\rightarrow q \bar{q}$"]
y_axis_BB_output = [0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725,
                    0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975
             ]
x_axis_Continuum_output = [0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725,
                    0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975
            ]


df_cm = pd.DataFrame(confusion_matrix.values, index = [i for i in y_axis_BB_output], columns = [i for i in x_axis_Continuum_output])
#df_cm = pd.DataFrame(confusion_matrix)

plt.figure(figsize = (30,20))
ax = sn.heatmap(df_cm, annot=True, cmap='RdYlBu')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize = 18, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 18)
plt.savefig("Matrix_large.png", bbox_inches='tight')
