import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
from scipy.stats import norm
from scipy.optimize import fsolve
from scipy.stats import gaussian_kde
from scipy.stats import percentileofscore

def GetMedianUncertainty(mu_name, color):
    # Plot with Seaborn
    fig, ax = plt.subplots()
    sns.histplot(df[mu_name], bins=51, kde=True, ax=ax, color=color)

    # Step 1: Calculate KDE with scipy
    kde = gaussian_kde(df[mu_name], bw_method='scott')
    x_vals = np.linspace(df[mu_name].min(), df['mu1'].max(), 1000)
    kde_vals = kde(x_vals)

    # Calculate the uncertainties
    median_value = np.percentile(df[mu_name], 50)
    lower_bound = np.percentile(df[mu_name], 15.866)
    upper_bound = np.percentile(df[mu_name], 84.134)
    lower_uncertainty = median_value - lower_bound
    upper_uncertainty = upper_bound - median_value

    # Display results
    ax.axvline(median_value, color='red', linestyle='--', label=f'Median: {median_value:.2f}')
    ax.axvline(lower_bound, color='green', linestyle='--', label=f'Lower bound: {lower_bound:.2f}')
    ax.axvline(upper_bound, color='orange', linestyle='--', label=f'Upper bound: {upper_bound:.2f}')
    ax.fill_betweenx([0, max(kde_vals)], lower_bound, upper_bound, color='gray', alpha=0.3)
    ax.legend()

    plt.show()

    print(f"Median: {median_value:.2f}")
    print(f"Lower uncertainty (16th percentile): {lower_uncertainty:.2f}")
    print(f"Upper uncertainty (84th percentile): {upper_uncertainty:.2f}")

    return median_value, upper_uncertainty, lower_uncertainty

# =========================== calculate correlation matrix and added UL ===========================

# constant value
BR1 = 0.0000048514
BR2 = 0.0000085024
BR3 = 0.0000156653

# input values
mu1_central = 0.877543
mu1_positive_error = 2.506356
mu1_negative_error = 2.260658

mu2_central = 2.590420
mu2_positive_error = 4.780812
mu2_negative_error = 4.378192

mu3_central = 4.210590
mu3_positive_error = 11.367196
mu3_negative_error = 9.878447

# calculate error. Just average. Any good idea?
mu1_error = (mu1_positive_error + mu1_negative_error) / 2.0
mu2_error = (mu2_positive_error + mu2_negative_error) / 2.0
mu3_error = (mu3_positive_error + mu3_negative_error) / 2.0

# Load the CSV file
df = pd.read_csv('mu_fit_values.csv')

# Set up the figure layout for subplots
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle(f"1D and 2D Plots for $\mu$")

# Draw 1D plot
sns.histplot(df['mu1'], bins=51, kde=True, ax=axs[0, 0], color="skyblue")
axs[0, 0].set_title(f"1D Plot of $\mu_{1}$")

sns.histplot(df['mu2'], bins=51, kde=True, ax=axs[0, 1], color="salmon")
axs[0, 1].set_title(f"1D Plot of $\mu_{2}$")

sns.histplot(df['mu3'], bins=51, kde=True, ax=axs[0, 2], color="lightgreen")
axs[0, 2].set_title(f"1D Plot of $\mu_{3}$")

# 2D scatter plots for each variable pair
sns.scatterplot(x='mu1', y='mu2', data=df, ax=axs[1, 0], color="blue")
axs[1, 0].set_title(f"2D Plot of $\mu_{1}$ vs $\mu_{2}$")

sns.scatterplot(x='mu2', y='mu3', data=df, ax=axs[1, 1], color="red")
axs[1, 1].set_title(f"2D Plot of $\mu_{2}$ vs $\mu_{3}$")

sns.scatterplot(x='mu1', y='mu3', data=df, ax=axs[1, 2], color="green")
axs[1, 2].set_title(f"2D Plot of $\mu_{1}$ vs $\mu_{3}$")

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# Calculate the correlation matrix
correlation_matrix = df.corr()
print("Correlation Matrix:\n", correlation_matrix)

# get sum of BR
BR_sum = (BR1 * mu1_central) + (BR2 * mu2_central) + (BR3 * mu3_central)
BR_error = math.sqrt((BR1 * mu1_error) * (BR1 * mu1_error) + \
           (BR2 * mu2_error) * (BR2 * mu2_error) + \
           (BR3 * mu3_error) * (BR3 * mu3_error) + \
           2 * correlation_matrix.loc["mu1", "mu2"] * (BR1 * mu1_error) * (BR2 * mu2_error) + \
           2 * correlation_matrix.loc["mu1", "mu3"] * (BR1 * mu1_error) * (BR3 * mu3_error) + \
           2 * correlation_matrix.loc["mu2", "mu3"] * (BR2 * mu2_error) * (BR3 * mu3_error))

print("sum of BR = %.8f +- %.8f" % (BR_sum, BR_error))

# Calculate UL
integral_from_minus_infty_to_zero = norm.cdf(0, loc=BR_sum, scale=BR_error)
target_integral = 0.9 * (1 - integral_from_minus_infty_to_zero)

# Define the function to solve for `u`
def equation(u):
    return norm.cdf(u, loc=BR_sum, scale=BR_error) - integral_from_minus_infty_to_zero - target_integral

BR_UL = fsolve(equation, x0=BR_sum)  # Start the search around the mean
print("UL of BR = %.6f" % BR_UL)

# =================================================================================================


# =========================== add UL directly from Toy MC distribution ============================

mu1_central_TOY, mu1_positive_error_TOY, mu1_negative_error_TOY = GetMedianUncertainty("mu1", "skyblue")
mu2_central_TOY, mu2_positive_error_TOY, mu2_negative_error_TOY = GetMedianUncertainty("mu2", "salmon")
mu3_central_TOY, mu3_positive_error_TOY, mu3_negative_error_TOY = GetMedianUncertainty("mu3", "lightgreen")

df['mu1_correct'] = (df['mu1'] - mu1_central_TOY) * (mu1_positive_error + mu1_negative_error) / (mu1_positive_error_TOY + mu1_negative_error_TOY) + mu1_central
df['mu2_correct'] = (df['mu2'] - mu2_central_TOY) * (mu2_positive_error + mu2_negative_error) / (mu2_positive_error_TOY + mu2_negative_error_TOY) + mu2_central
df['mu3_correct'] = (df['mu3'] - mu3_central_TOY) * (mu3_positive_error + mu3_negative_error) / (mu3_positive_error_TOY + mu3_negative_error_TOY) + mu3_central

df['BR'] = (df['mu1_correct'] * BR1) + (df['mu2_correct'] * BR2) + (df['mu3_correct'] * BR3)

GetMedianUncertainty("BR", "indigo")

percentile_rank = percentileofscore(df["BR"], 0)
upper_limit = np.percentile(df["BR"], percentile_rank + (100 - percentile_rank) * 0.9)
print("UL of BR = %.6f" % upper_limit)

# =================================================================================================

# ========================================== sanity check =========================================

mu1_central_TOY, mu1_positive_error_TOY, mu1_negative_error_TOY = GetMedianUncertainty("mu1_correct", "skyblue")
mu2_central_TOY, mu2_positive_error_TOY, mu2_negative_error_TOY = GetMedianUncertainty("mu2_correct", "salmon")
mu3_central_TOY, mu3_positive_error_TOY, mu3_negative_error_TOY = GetMedianUncertainty("mu3_correct", "lightgreen")

# =================================================================================================
