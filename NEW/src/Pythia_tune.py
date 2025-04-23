import numpy as np
import pandas as pd
import os
from glob import glob

# Reference data and uncertainties (3x9 arrays)
ref_Xs_gamma = 0.01 * np.array([
    [10.6, 5.88, 23.2, 44.5, 0.46, 9.99, 0.52, 4.78, 0.00],
    [3.12, 1.13, 15.7, 20.6, 9.48, 26.9, 5.29, 15.5, 2.20],
    [1.65, 0.51, 9.02, 7.98, 16.8, 26.7, 10.3, 21.8, 5.38]
])

rel_unc_Xs_gamma = np.array([
    [0.0651, 0.0990, 0.0554, 0.0400, 0.6304, 0.1175, 0.7115, 0.2395, 0.0000],
    [0.1053, 0.2727, 0.0686, 0.0741, 0.0942, 0.0733, 0.2777, 0.1715, 0.1168],
    [0.4250, 1.4490, 0.2646, 0.4003, 0.2509, 0.2671, 0.8409, 0.6682, 0.2898]
])

abs_unc_Xs_gamma = rel_unc_Xs_gamma * ref_Xs_gamma

ref_Xs_gamma_norm = ref_Xs_gamma[:2, :7] / ref_Xs_gamma[:2, :7].sum(axis=1, keepdims=True)

ref_Xs_Jpsi = np.array([
    [0.13, 0.07, 0.30, 0.30, 0.05, 0.15, 0.01, 0.00],
    [0.06, 0.03, 0.18, 0.19, 0.13, 0.27, 0.12, 0.02]
])

abs_unc_Xs_Jpsi = np.array([
    [0.01, 0.01, 0.03, 0.04, 0.06, 0.09, 0.11, 0.01],
    [0.01, 0.01, 0.03, 0.04, 0.06, 0.09, 0.11, 0.01]
])

ref_Xs_Jpsi_norm = ref_Xs_Jpsi[:2, :7] / ref_Xs_Jpsi[:2, :7].sum(axis=1, keepdims=True)


# Load CSVs and compute chi² from the uncertainty of Xs gamma study only
def compute_chi2_and_sigmas(file):
    df = pd.read_csv(file, header=None)
    data = df.iloc[:3, :9].values
    diff = np.abs(data - ref_Xs_gamma)
    mask = abs_unc_Xs_gamma != 0  # Ignore zero-uncertainty entries

    chi2 = np.sum((diff ** 2)[mask] / (abs_unc_Xs_gamma ** 2)[mask])

    within_1sigma = np.sum((diff <= 1 * abs_unc_Xs_gamma)[mask])
    within_2sigma = np.sum((diff <= 2 * abs_unc_Xs_gamma)[mask])
    within_3sigma = np.sum((diff <= 3 * abs_unc_Xs_gamma)[mask])
    total_valid = np.sum(mask)

    return chi2, within_1sigma, within_2sigma, within_3sigma, total_valid

# Load CSVs and compute chi² . Here, sigma is the difference between Xsgamma and XsJpsi
def compute_chi2_and_sigmas_from_Xsgamma(file):
    df = pd.read_csv(file, header=None)
    data = df.iloc[:2, :7].values
    data_norm = data / data[:2, :7].sum(axis=1, keepdims=True)
    diff = np.abs(data_norm - ref_Xs_gamma_norm)
    mask = abs(ref_Xs_gamma_norm - ref_Xs_Jpsi_norm) != 0  # Ignore zero-uncertainty entries

    chi2 = np.sum((diff ** 2)[mask] / ((ref_Xs_gamma_norm - ref_Xs_Jpsi_norm) ** 2)[mask])

    within_1sigma = np.sum((diff <= 1 * abs(ref_Xs_gamma_norm - ref_Xs_Jpsi_norm))[mask])
    within_2sigma = np.sum((diff <= 2 * abs(ref_Xs_gamma_norm - ref_Xs_Jpsi_norm))[mask])
    within_3sigma = np.sum((diff <= 3 * abs(ref_Xs_gamma_norm - ref_Xs_Jpsi_norm))[mask])
    total_valid = np.sum(mask)

    return chi2, within_1sigma, within_2sigma, within_3sigma, total_valid

# Directory with CSVs
csv_dir = "./"  # Change this if needed
csv_files = sorted(glob(os.path.join(csv_dir, "*.csv")))

# Compute and sort by chi2
results = [(f, *compute_chi2_and_sigmas_from_Xsgamma(f)) for f in csv_files]
results.sort(key=lambda x: x[1])  # sort by chi²

for fname, chi2, w1, w2, w3, total in results:
    print(f"{os.path.basename(fname)}: chi² = {chi2:.4f}, "
          f"within 1σ = {w1}/{total}, 2σ = {w2}/{total}, 3σ = {w3}/{total}")
