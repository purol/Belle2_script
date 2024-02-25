import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation

# Define a function to plot and save a dataset
def plot_and_save_data(x, y1, yerr1, filename, y1_label, title):
    fig, ax = plt.subplots(figsize=(8, 7))

    trans1 = ax.transData + ScaledTranslation(-5/72, 0, fig.dpi_scale_trans)
    er1 = ax.errorbar(x, y1, yerr=yerr1, marker="o", linestyle="none", transform=trans1, label=y1_label)

    # Rotate the x-axis labels by 45 degrees
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=30, ha='right')

    # Add a title and labels for the axes
    ax.set_xlabel('Decay Mode')
    ax.set_ylabel('FBDT efficiency')

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

x = [r'$K$', r'$K \pi$', r'$K 2\pi$', r'$K 3\pi$', r'$K 4\pi$', r'including $3K$']

y1_1 = [0.6174882586281491, 0.0872938566356157, 0.0110777153357214, 0, 0, 0]
yerr1_1 = [0.09172157909 * 0.6174882586281491, 0.19349504716 * 0.0872938566356157, 1.1099267076 * 0.0110777153357214, 0, 0, 0]

plot_and_save_data(x, y1_1, yerr1_1, "FBDT_eff.png", r'$B \rightarrow X_{s} J/\psi$', r'')
