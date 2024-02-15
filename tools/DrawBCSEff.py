import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation

# Define a function to plot and save a dataset
def plot_and_save_data(x, y1, yerr1, y2, yerr2, filename, y1_label, y2_label, title):
    fig, ax = plt.subplots(figsize=(8, 7))

    trans1 = ax.transData + ScaledTranslation(-5/72, 0, fig.dpi_scale_trans)
    trans2 = ax.transData + ScaledTranslation(+5/72, 0, fig.dpi_scale_trans)
    er1 = ax.errorbar(x, y1, yerr=yerr1, marker="o", linestyle="none", transform=trans1, label=y1_label)
    er2 = ax.errorbar(x, y2, yerr=yerr2, marker="o", linestyle="none", transform=trans2, label=y2_label)

    # Rotate the x-axis labels by 45 degrees
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=30, ha='right')

    # Add a title and labels for the axes
    ax.set_xlabel('Decay Mode')
    ax.set_ylabel('Fraction of decay modes[%]')

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

x = [r'$K$', r'$K\pi$', r'$K2\pi$', r'$K3\pi$', r'$K4\pi$', r'$3K$', r'$3K\pi$']

y1_1 = [39.951935, 32.761871, 12.137244, 1.756578, 0.15261, 0.185272, 0.014101]
total_sum = sum(y1_1)
y1_1 = [(value / total_sum) * 100 for value in y1_1]
yerr1_1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

y2_1 = [24.506693, 19.15963, 6.137504, 0.801534, 0.046628, 0.109851, 0.008537]
total_sum = sum(y2_1)
y2_1 = [(value / total_sum) * 100 for value in y2_1]
yerr2_1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

plot_and_save_data(x, y1_1, yerr1_1, y2_1, yerr2_1, "fraction_modes.png", r'Before BCS', r'After BCS', r'')
