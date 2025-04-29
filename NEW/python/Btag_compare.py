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
    ax.set_ylim(-1.0, y_max * 1.2)
    ax.set_xlim(-0.5, len(x) - 0.5)

    # Rotate the x-axis labels by 45 degrees
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=30, ha='right', fontsize=14)

    # Add a title and labels for the axes
    ax.set_xlabel(r'$B_{tag}$ Decay Mode', fontsize=14)
    ax.set_ylabel('The number of event', fontsize=14)

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

y1 = [10.793207, 19.589011, 24.772498, 23.985855, 7.603868, 6.379483, 9.534895, 9.217291, 2.823215, 7.177456, 4.093737, 16.610551] # for MC
yerr1 = [1.592359, 2.152339, 2.419511, 2.377383, 1.345935, 1.228819, 1.509222, 1.476520, 0.815984, 1.311465, 0.994334, 1.986981]

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

y1 = [1.162826, 2.836806, 5.274539, 3.647890, 1.641918, 1.654976, 4.735004, 4.100285, 5.223389, 0.244711, 2.869994] # for MC
yerr1 = [0.520266, 0.819345, 1.125479, 0.946237, 0.620797, 0.628498, 1.059227, 0.995556, 1.115413, 0.244711, 0.828931]

y2 = [2, 4, 3, 8, 1, 3, 0, 4, 2, 2, 3] # for data
yerr2 = [1.414214, 2.000000, 1.732051, 2.828427, 1.000000, 1.732051, 0, 2.000000, 1.414214, 1.414214, 1.732051]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "Number_of_event_comparison_Bzero", r'$B \rightarrow K^{(*)} J/\psi$ MC', r'$B \rightarrow K^{(*)} J/\psi$ data', r'')
