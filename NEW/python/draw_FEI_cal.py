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
e_values = [0.8562, 1.1172, 0.8961, 0.8518, 0.7957, 1.0491, 1.0418, 0.8610, 0.7607, 1.0084, 0.7157, 0.4609, 0.8197]
e_errors = [0.0321, 0.0431, 0.0344, 0.0323, 0.0303, 0.0409, 0.0409, 0.0335, 0.0298, 0.0417, 0.0289, 0.0193, 0.0314]

mu_values = [0.8276, 1.0691, 0.8640, 0.8297, 0.7671, 1.0288, 1.0424, 0.8819, 0.7204, 1.0276, 0.7023, 0.4668, 0.8046]
mu_errors = [0.0311, 0.0417, 0.0335, 0.0317, 0.0295, 0.0405, 0.0412, 0.0345, 0.0286, 0.0434, 0.0292, 0.0199, 0.0309]

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
e_values = [0.8956, 1.1425, 1.0372, 0.9459, 0.9285, 1.0833, 1.0343, 0.7022, 0.7326, 0.8008, 0.6215, 0.9371]
e_errors = [0.0330, 0.0468, 0.0406, 0.0362, 0.0356, 0.0426, 0.0428, 0.0276, 0.0284, 0.0308, 0.0290, 0.0358]

mu_values = [0.8796, 1.1661, 1.0303, 0.9179, 0.8941, 1.0960, 1.0569, 0.6922, 0.7169, 0.7825, 0.6104, 0.9209]
mu_errors = [0.0325, 0.0479, 0.0410, 0.0358, 0.0348, 0.0436, 0.0442, 0.0273, 0.0281, 0.0304, 0.0298, 0.0354]

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
