#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import FixedLocator

cov_Bp = np.load('FEI_cal_covariance_Bp_001.npy')
cov_B0 = np.load('FEI_cal_covariance_B0_001.npy')

np.savetxt("FEI_cal_Bp_cov_all.txt",cov_Bp,fmt='%1.12f')
np.savetxt("FEI_cal_B0_cov_all.txt",cov_B0,fmt='%1.12f')
