#!/usr/bin/env python3

import os
import sys
import json
import warnings
import itertools
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
import math

plot_title = 'e$^{-}$ ratio table for cut "kaonID>0.6"'
fig, ax = plt.subplots(1, 1, figsize=(7, 5), dpi=120)

df = pd.read_csv("./kaonID_efficiency_electron.csv")

x_label_name = "theta"
y_label_name = "p"
x_label = [0.22, 0.56, 1.13, 1.57, 1.88, 2.23, 2.71] # theta
y_label = [0.4, 0.5, 1.0, 1.5, 2, 2.5, 3] # p

arr_charge = df["charge"].to_numpy()
arr_p_min = df["p_min"].to_numpy()
arr_p_max = df["p_max"].to_numpy()
arr_theta_min = df["theta_min"].to_numpy()
arr_theta_max = df["theta_max"].to_numpy()
arr_data_MC_ratio = df["data_MC_ratio"].to_numpy()
arr_data_MC_uncertainty_stat_up = df["data_MC_uncertainty_stat_up"].to_numpy()
arr_data_MC_uncertainty_stat_dn = df["data_MC_uncertainty_stat_dn"].to_numpy()
arr_data_MC_uncertainty_sys_up = df["data_MC_uncertainty_sys_up"].to_numpy()
arr_data_MC_uncertainty_sys_dn = df["data_MC_uncertainty_sys_dn"].to_numpy()
arr_data_MC_uncertainty_up = np.sqrt(arr_data_MC_uncertainty_stat_up*arr_data_MC_uncertainty_stat_up+arr_data_MC_uncertainty_sys_up*arr_data_MC_uncertainty_sys_up)
arr_data_MC_uncertainty_dn = np.sqrt(arr_data_MC_uncertainty_stat_dn*arr_data_MC_uncertainty_stat_dn+arr_data_MC_uncertainty_sys_dn*arr_data_MC_uncertainty_sys_dn)

arr = np.empty((len(x_label) - 1,len(y_label) - 1,))
arr[:] = np.nan

annot_labels = [["" for col in range(len(x_label) - 1)] for row in range(len(y_label) - 1)]

for charge, p_min, p_max, theta_min, theta_max, data_MC_ratio, data_MC_uncertainty_up, data_MC_uncertainty_dn in zip(arr_charge, arr_p_min, arr_p_max, arr_theta_min, arr_theta_max, arr_data_MC_ratio, arr_data_MC_uncertainty_up, arr_data_MC_uncertainty_dn):
    if charge == "+":
        continue

    p_index = 0
    theta_index = 0
    for i, temp_p_min in zip(range(0,len(y_label)), y_label):
        if math.isclose(temp_p_min, p_min):
            p_index = i
            break
    for i, temp_theta_min in zip(range(0,len(x_label)), x_label):
        if math.isclose(temp_theta_min, theta_min):
            theta_index = i
            break
    arr[p_index][theta_index] = data_MC_ratio
    annot_labels[p_index][theta_index] = f'{data_MC_ratio:.2f}$^{{+{data_MC_uncertainty_up:.2f}}}_{{-{data_MC_uncertainty_dn:.2f}}}$'


plot_args = {'linewidths': .5, 'fmt': '', 'cbar_kws': {'fraction': .05}, 'vmin': 0, 'vmax': 1, 'cmap': 'plasma'}


sns.heatmap(arr, annot=annot_labels, ax=ax, **plot_args)

ax.set_xticks(np.arange(0, len(x_label)))
ax.set_xticklabels(np.round(x_label, 4))
ax.set_yticks(np.arange(0, len(y_label)))
ax.set_yticklabels(np.round(y_label, 4), rotation=0)
ax.set_xlabel(f'{x_label_name} bins', fontsize=14)
ax.set_ylabel(f'{y_label_name} bins', fontsize=14)

ax.set_title(plot_title)

plt.savefig('E_minus_kaonID.png')
