import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation

# Define a function to plot and save a dataset
def plot_and_save_data(x, y1, yerr1, y2, yerr2, filename, y1_label, y2_label, title):
    fig, ax = plt.subplots(figsize=(8, 7))

    trans1 = ax.transData + ScaledTranslation(-5/72, 0, fig.dpi_scale_trans)
    trans2 = ax.transData + ScaledTranslation(+5/72, 0, fig.dpi_scale_trans)
    er1 = ax.errorbar(x, y1, yerr=yerr1, marker="o", linestyle="none", transform=trans1, label=y1_label)
    er2 = ax.errorbar(x, y2, yerr=yerr2, marker="v", linestyle="none", transform=trans2, label=y2_label)

    # Rotate the x-axis labels by 45 degrees
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=30, ha='right')

    # Add a title and labels for the axes
    ax.set_xlabel('Decay Mode')
    ax.set_ylabel('calibration factor')

    ax.legend(loc='best')

    ax.text(0.02, 0.95, title, transform=ax.transAxes, color='black', fontsize=16)

    plt.tight_layout()

    # Save the figure as a PNG image
    plt.savefig(filename, dpi=300)
    #plt.show()
    plt.close()  # Close the figure to release resources

def print_ratio(y1, yerr1, y2, filename):
    ratios = [y1_value / y2_value for y1_value, y2_value in zip(y1, y2)]
    relative_uncertainties = [yerr1_value / y1_value if y1_value != 0 else 0 for yerr1_value, y1_value in zip(yerr1, y1)]
    with open(filename, "w") as file:
        file.write("%d\n"%len(y1))
        for y1_value, y2_value, ratio, relative_uncertainty in zip(y1, y2, ratios, relative_uncertainties):
            file.write(f"{y1_value} {y2_value} {ratio} {relative_uncertainty}\n")

x = [r'$\bar{D}^{0}\pi^{+}$', \
     r'$\bar{D}^{0}\pi^{+}\pi^{0}$', \
     r'$\bar{D}^{0}\pi^{+}\pi^{+}\pi^{-}$', \
     r'$\bar{D}^{0}\pi^{+}\pi^{+}\pi^{-}\pi^{0}$', \
     r'$\bar{D}^{0*}\pi^{+}$', \
     r'$\bar{D}^{0*}\pi^{+}\pi^{0}$', \
     r'$\bar{D}^{0*}\pi^{+}\pi^{+}\pi^{-}$', \
     r'$\bar{D}^{0*}\pi^{+}\pi^{+}\pi^{-}\pi^{0}$', \
     r'$D^{-}\pi^{+}\pi^{+}$', \
     r'$D^{-}\pi^{+}\pi^{+}\pi^{0}$', \
     r'$\Lambda_{c}^{-}p\pi^{+}\pi^{-}\pi^{+}$', \
     'rest']

y1 = [0.956417, 0.867980, 0.620128, 0.724721, 0.770419, 0.370744, 0.880669, 0.570453, 0.615938, 0.225165, 0.0, 0.073675] # for B->K J/psi
yerr1 = [0.373878, 0.254466, 0.191536, 0.226227, 0.418873, 0.273338, 0.395046, 0.303658, 0.658487, 0.231018, 0.0, 0.074309]

y2 = [1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59] # for BELLE2-NOTE-PH-2023-008
yerr2 = [0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "FEI_calibration_compare", r'$B \rightarrow K J/\psi$ study', 'BELLE2-NOTE-PH-2023-008', r'')
