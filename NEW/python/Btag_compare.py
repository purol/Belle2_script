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

y1 = [10.767046, 19.541529, 24.712453, 23.927716, 7.585437, 6.364019, 9.511783, 9.194950, 2.816372, 7.160058, 4.083815, 16.570289] # for MC
yerr1 = [1.588499, 2.147122, 2.413646, 2.371621, 1.342673, 1.225840, 1.505564, 1.472941, 0.814006, 1.308287, 0.991924, 1.982165]

y2 = [14, 18, 19, 18, 5, 4, 6, 4, 1, 3, 0, 3] # for data
yerr2 = [3.741657, 4.242641, 4.358899, 4.242641, 2.236068, 2.000000, 2.449490, 2.000000, 1.000000, 1.732051, 0.000000, 1.732051]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "Number_of_event_comparison_Bplus", r'$B \rightarrow K^{(*)} J/\psi$ MC', r'$B \rightarrow K^{(*)} J/\psi$ data', r'')

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

y1 = [1.160007, 2.829930, 5.261755, 3.639048, 1.637938, 1.650964, 4.723527, 4.090346, 5.210728, 0.244118, 2.863038] # for MC
yerr1 = [0.519005, 0.817360, 1.122751, 0.943944, 0.619292, 0.626974, 1.056660, 0.993143, 1.112710, 0.244118, 0.826922]

y2 = [2, 4, 3, 8, 1, 3, 0, 4, 2, 2, 3] # for data
yerr2 = [1.414214, 2.000000, 1.732051, 2.828427, 1.000000, 1.732051, 0, 2.000000, 1.414214, 1.414214, 1.732051]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "Number_of_event_comparison_Bzero", r'$B \rightarrow K^{(*)} J/\psi$ MC', r'$B \rightarrow K^{(*)} J/\psi$ data', r'')
