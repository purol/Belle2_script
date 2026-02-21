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

    y1_top = [y + err for y, err in zip(y1, yerr1)]
    y2_top = [y + err for y, err in zip(y2, yerr2)]
    y_max = max(max(y1_top), max(y2_top))
    ax.set_ylim(0.0, y_max * 1.2)
    ax.set_xlim(-0.5, len(x) - 0.5)

    # Rotate the x-axis labels by 45 degrees
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=30, ha='right', fontsize=14)

    # Add a title and labels for the axes
    ax.set_xlabel(r'$B_{tag}$ Decay Mode', fontsize=14)
    ax.set_ylabel('FEI calibration factor', fontsize=14)

    ax.legend(loc='best', fontsize=14)

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

# for Bplus
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

y1 = [1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59] # for MC15ri
yerr1 = [0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02]

y2 = [1.0563843143394622, 0.824350934769151, 0.7821931102725646, 0.6720352787461095, 1.003394557399459, 1.0056643773735399, 0.8369671783480447, 0.6783049826281483, 0.88198076213399, 0.6437129465475984, 0.3840523705585905, 0.7727764587488754] # for MC15rd
yerr2 = [0.029112659708673237, 0.061272260861712384, 0.0686254854670932, 0.1594569267928963, 0.028165156667679287, 0.028379916542260016, 0.027153130014423914, 0.08008729127106695, 0.1985684611368228, 0.09017947779294383, 0.12049932626953824, 0.03556587244419456]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "charged_B_FEI_cal", r'MC15ri', r'MC15rd', r'')

# for B zero
x = [r'$D^{-}\pi^{+}$', \
     r'$D^{-}\pi^{+}\pi^{0}$', \
     r'$D^{-}\pi^{+}\pi^{+}\pi^{-}$', \
     r'$D^{-}\pi^{+}\pi^{+}\pi^{-}\pi^{0}$', \
     r'$\bar{D}^{0}\pi^{+}\pi^{-}$', \
     r'$D^{-*}\pi^{+}$', \
     r'$D^{-*}\pi^{+}\pi^{0}$', \
     r'$D^{-*}\pi^{+}\pi^{+}\pi^{-}$', \
     r'$D^{-*}\pi^{+}\pi^{+}\pi^{-}\pi^{0}$', \
     r'$\Lambda_{c}^{-}p\pi^{+}\pi^{-}$', \
     'rest']

y1 = [1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79] # for MC15ri
yerr1 = [0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12]

y2 = [1.1351815678693125, 1.0152536119073212, 0.900827789502242, 0.820330582165554, 0.9601307310429825, 1.0287495457230427, 0.6565638369059716, 0.6884022887870119, 0.7376123919921236, 0.6379774334598569, 0.9912539950874157] # for MC15rd
yerr2 = [0.03616790985932927, 0.031621512590247744, 0.027763030442755705, 0.026388965159170368, 0.03181259001191749, 0.03308978800268333, 0.020328451739706932, 0.021159300223306747, 0.022951063339037594, 0.026623990893496787, 0.08285314589226657]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "neutral_B_FEI_cal", r'MC15ri', r'MC15rd', r'')
