#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt


confusion_matrix = pd.read_excel("./FOM_small.xlsx", header=None, index_col=None, nrows=20)
#print(confusion_matrix)
confusion_matrix.to_numpy()

#axis_list = [r"$\tau^{\mp}\rightarrow\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{\pm}\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\nu$",r"$\tau^{\mp}\rightarrow\pi^{\mp}\pi^{0}\pi^{0}\nu$",r"$Z\rightarrow q \bar{q}$"]
y_axis_BB_output = [0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
                    0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49
             ]
x_axis_Continuum_output = [0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
                    0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49
            ]


df_cm = pd.DataFrame(confusion_matrix.values, index = [i for i in y_axis_BB_output], columns = [i for i in x_axis_Continuum_output])
#df_cm = pd.DataFrame(confusion_matrix)

plt.figure(figsize = (30,20))
ax = sn.heatmap(df_cm, annot=True, cmap='RdYlBu')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize = 18, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize = 18)
plt.savefig("Matrix_small.png", bbox_inches='tight')
