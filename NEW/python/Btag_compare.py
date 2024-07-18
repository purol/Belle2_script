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
    ax.set_ylabel('The number of event')

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

y1 = [8.364552, 16.129407, 19.350836, 16.558090, 5.191980, 5.394551, 6.813006, 7.011976, 1.623541, 4.441195, 3.351034, 13.573115] # for MC
yerr1 = [1.395016, 1.943602, 2.125431, 1.966790, 1.108748, 1.125932, 1.266428, 1.280623, 0.613804, 1.019269, 0.897141, 1.783717]

y2 = [8, 14, 12, 12, 4, 2, 6, 4, 1, 1, 0, 1] # for data
yerr2 = [2.828427, 3.741657, 3.464102, 3.464102, 2.0, 1.414214, 2.449490, 2, 1, 1, 0, 1]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "Number_of_event_comparison", r'$B \rightarrow K J/\psi$ MC', r'$B \rightarrow K J/\psi$ data', r'')
