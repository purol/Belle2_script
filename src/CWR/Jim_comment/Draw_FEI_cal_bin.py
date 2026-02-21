import numpy as np
import matplotlib.pyplot as plt

# Data
SIGNAL_Nevt = np.array([1.686468, 4.500264, 1.385178, 2.185774, 4.60187,
                        0.603308, 0.831906, 0.890514, 1.464753, 2.732566,
                        0.703914, 1.033782, 0.761144, 1.184804, 1.585798])
BKG_Nevt = np.array([570.6563, 638.5276, 81.35691, 73.31551, 44.94842,
                     114.6054, 111.087, 81.48235, 78.17103, 46.28284,
                     509.0142, 486.8386, 220.2566, 195.4779, 85.98803])

SIGNAL_Nevt_without_FEI = np.array([2.136508, 5.668246, 1.738668, 2.738419, 5.626477,
                                    0.749239, 1.023037, 1.092724, 1.813467, 3.280265,
                                    0.874265, 1.283692, 0.94099, 1.455113, 1.922771])
BKG_Nevt_without_FEI = np.array([620.5895, 709.9028, 91.30438, 84.35691, 51.57401,
                                  127.8642, 124.3255, 91.52758, 90.19963, 54.27302,
                                  587.2955, 566.031, 257.6299, 230.0471, 100.2072])

# x-axis indices
bins = np.arange(1, len(SIGNAL_Nevt) + 1)

# Compute ratios
signal_ratio = SIGNAL_Nevt / SIGNAL_Nevt_without_FEI
bkg_ratio = BKG_Nevt / BKG_Nevt_without_FEI

print("Signal ratio (with/without FEI):")
print(signal_ratio)
print("\nBackground ratio (with/without FEI):")
print(bkg_ratio)

def make_ratio_plot(x, y1, y2, ylabel, title):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, 
                                   gridspec_kw={'height_ratios': [3, 1]},
                                   figsize=(8, 6))
    # Top plot: two lines
    ax1.bar(x, y2, width=1.0, facecolor='none', hatch='////', label='Without FEI calibration', edgecolor='red')
    ax1.bar(x, y1, width=1.0, facecolor='none', hatch='\\\\', label='With all calibration', edgecolor='blue')
    ax1.set_ylabel(ylabel, fontsize=14)
    ax1.legend(fontsize=12)
    ax1.set_title(title, fontsize=16)

    ax1.set_xticks(x)
    ax2.set_xticks(x)

    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)

    # Bottom plot: ratio
    ratio = y1 / y2
    ax2.plot(x, ratio, 'ko')
    ax2.axhline(1.0, color='red', linestyle='--')
    ax2.set_ylabel('Ratio', fontsize=14)
    ax2.set_xlabel('Bin index', fontsize=14)
    ax2.set_ylim(0.6, 1.1)

    plt.tight_layout()
    plt.show()

# Plot for signal
make_ratio_plot(bins, SIGNAL_Nevt, SIGNAL_Nevt_without_FEI,
                ylabel='Events', title='Signal Events')

# Plot for background
make_ratio_plot(bins, BKG_Nevt, BKG_Nevt_without_FEI,
                ylabel='Events', title='Background Events')
