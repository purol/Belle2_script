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

y1 = [10.684368, 19.391474, 24.522691, 23.743979, 7.527190, 6.315151, 9.438744, 9.124343, 2.794746, 7.105078, 4.052456, 16.443049] # for MC
yerr1 = [1.576302, 2.130635, 2.395113, 2.353409, 1.332362, 1.216427, 1.494003, 1.461631, 0.807755, 1.298241, 0.984307, 1.966945]

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

y1 = [1.151100, 2.808200, 5.221351, 3.611105, 1.625361, 1.638287, 4.687256, 4.058938, 5.170716, 0.242243, 2.841053] # for MC
yerr1 = [0.515020, 0.811083, 1.114129, 0.936695, 0.614536, 0.622160, 1.048546, 0.985517, 1.104165, 0.242243, 0.820572]

y2 = [2, 4, 3, 8, 1, 3, 0, 4, 2, 2, 3] # for data
yerr2 = [1.414214, 2.000000, 1.732051, 2.828427, 1.000000, 1.732051, 0, 2.000000, 1.414214, 1.414214, 1.732051]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "Number_of_event_comparison_Bzero", r'$B \rightarrow K^{(*)} J/\psi$ MC', r'$B \rightarrow K^{(*)} J/\psi$ data', r'')
