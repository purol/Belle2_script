#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

y = np.arange(35)
decay = [r"$B^{+}\rightarrow K^{+}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{*+}\rightarrow K^{+}\pi^{0}\nu\bar{\nu}$",
            r"$B^{+}\rightarrow K^{*+}\rightarrow K^{0}\pi^{+}\nu\bar{\nu}$",      
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
            r"$B^{0}\rightarrow K^{*0}\rightarrow K^{+}\pi^{-}\nu\bar{\nu}$",
            r"$B^{0}\rightarrow K^{*0}\rightarrow K^{0}\pi^{0}\nu\bar{\nu}$",
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
            r"$others$"]

Nevt = [7.194184,
        4.059653,
        0.981578,
        1.733003,
        0.606261,
        0.605919,
        0.256311,
        0.268614,
        0.005810,
        0.001025,
        0.005468,
        0.287068,
        0.013328,
        0.026998,
        0.058439,
        0.002392,
        0.007860,
        0.650399,
        2.624511,
        0.385177,
        1.385979,
        0.207845,
        0.673859,
        0.089963,
        0.018613,
        0.027575,
        0.018613,
        0.000345,
        0.030332,
        0.053771,
        0.004136,
        0.005515,
        0.000689,
        0.000000,
        0.662485]

decay.reverse()
Nevt.reverse()

plt.rcParams['figure.figsize'] = [15, 12]

fig, ax = plt.subplots()
bars = ax.barh(decay, Nevt)
#ax.bar_label(bars)

#plt.barh(y, Nevt)
#plt.yticks(y, decay)

plt.xlabel("Number of Event (MC level)")
plt.ylabel("Decay Modes (MC level)")
plt.grid(True, axis = "x", linestyle="--")

plt.savefig("Nevt.png", bbox_inches='tight')
