import matplotlib.pyplot as plt
import numpy as np

# ========================= for B plus ========================= 
labels = [
    r"All", r"$\bar{D}^0\pi^+$", r"$\bar{D}^0\pi^+\pi^0$", 
    r"$\bar{D}^0\pi^+\pi^+\pi^-$", r"$\bar{D}^0\pi^+\pi^+\pi^-\pi^0$", 
    r"$\bar{D}^{*0}\pi^+$", r"$\bar{D}^{*0}\pi^+\pi^0$", 
    r"$\bar{D}^{*0}\pi^+\pi^+\pi^-$", r"$\bar{D}^{*0}\pi^+\pi^+\pi^-\pi^0$", 
    r"$D^-\pi^+\pi^+$", r"$D^-\pi^+\pi^+\pi^0$", 
    r"$\Lambda_c^{-} p\pi^+\pi^-\pi^+$", "rest"
]
e_values = [0.734, 1.083, 0.8368, 0.7488, 0.6887, 1.009, 0.9859, 0.7596, 0.6411, 0.9597, 0.6471, 0.3934, 0.6396]
e_errors = [0.0213, 0.0344, 0.0258, 0.0231, 0.0209, 0.0332, 0.0316, 0.0239, 0.0203, 0.036, 0.022, 0.0145, 0.0188]

mu_values = [0.7205, 1.067, 0.8111, 0.7326, 0.6499, 1.001, 0.9809, 0.7734, 0.6009, 1.001, 0.6423, 0.4117, 0.5958]
mu_errors = [0.0206, 0.0336, 0.0247, 0.0221, 0.0196, 0.0323, 0.0313, 0.024, 0.019, 0.0369, 0.0227, 0.0151, 0.0179]

# Reverse the labels and values for horizontal plot
labels = labels[::-1]
e_values = e_values[::-1]
e_errors = e_errors[::-1]
mu_values = mu_values[::-1]
mu_errors = mu_errors[::-1]

y_positions = np.arange(len(labels))

# Create the plot
plt.figure(figsize=(8, 6))

plt.errorbar(e_values, y_positions, xerr=e_errors, fmt='o', color='black', label=r"$e$")
plt.errorbar(mu_values, y_positions, xerr=mu_errors, fmt='o', color='blue', label=r"$\mu$")

plt.yticks(ticks=y_positions, labels=labels)
plt.xlabel("Calibration Factor")
plt.xlim(0, 1.4)
plt.legend(loc="lower right")
plt.title(r"$\mathcal{P} > 0.001, B^+$")
plt.tight_layout()

# Display the plot
plt.show()

# ========================= for B zero ========================= 
labels = [
    r"$B^0_{tag}e^-$", r"$D^-\pi^+$", r"$D^-\pi^+\pi^0$", 
    r"$D^-\pi^+\pi^+\pi^-$", r"$D^-\pi^+\pi^+\pi^-\pi^0$", 
    r"$\bar{D}^0\pi^+\pi^-$", r"$D^{*-}\pi^+$", 
    r"$\bar{D}^{*-}\pi^+\pi^0$", r"$\bar{D}^{*-}\pi^+\pi^+\pi^-$", 
    r"$\bar{D}^{*-}\pi^+\pi^+\pi^-\pi^0$", 
    r"$\Lambda_c^- p\pi^+\pi^-$", "rest"
]
e_values = [0.7956, 1.173, 0.9856, 0.8693, 0.7991, 0.9909, 1.029, 0.6749, 0.6691, 0.7115, 0.4898, 0.6993]
e_errors = [0.0233, 0.0451, 0.0338, 0.0286, 0.0266, 0.0348, 0.0393, 0.0223, 0.0215, 0.0226, 0.0234, 0.022]

mu_values = [0.7725, 1.176, 0.9413, 0.8354, 0.7882, 1.004, 1.013, 0.6718, 0.6493, 0.6825, 0.4918, 0.6815]
mu_errors = [0.0224, 0.0431, 0.032, 0.0274, 0.026, 0.0359, 0.039, 0.0222, 0.0211, 0.0217, 0.0227, 0.0212]

# Reverse the labels and values for horizontal plot
labels = labels[::-1]
e_values = e_values[::-1]
e_errors = e_errors[::-1]
mu_values = mu_values[::-1]
mu_errors = mu_errors[::-1]

y_positions = np.arange(len(labels))

# Create the plot
plt.figure(figsize=(8, 6))

plt.errorbar(e_values, y_positions, xerr=e_errors, fmt='o', color='black', label=r"$e$")
plt.errorbar(mu_values, y_positions, xerr=mu_errors, fmt='o', color='blue', label=r"$\mu$")

plt.yticks(ticks=y_positions, labels=labels)
plt.xlabel("Calibration Factor")
plt.xlim(0, 1.4)
plt.legend(loc="lower right")
plt.title(r"$\mathcal{P} > 0.001, B^0$")
plt.tight_layout()

# Display the plot
plt.show()
