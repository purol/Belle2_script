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
    ax.set_ylim(-0.3, y_max * 1.2)
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

x = ['category 1','category 2','category 3']

y1 = [78.342511, 32.405429, 30.395329] # for B->K(*) J/psi MC
yerr1 = [4.277783, 2.761032, 2.678736]

y2 = [69, 19, 7] # for data
yerr2 = [8.306624, 4.358899, 2.645751]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "Number_of_event_comparison_Bplus_My", r'$B \rightarrow K^{(*)} J/\psi$ MC', r'$B \rightarrow K^{(*)} J/\psi$ data', r'')

x = ['category 1']

y1 = [33.510046] # for B->K(*) J/psi MC
yerr1 = [2.836296]

y2 = [32] # for data
yerr2 = [5.656854]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "Number_of_event_comparison_Bzero_My", r'$B \rightarrow K^{(*)} J/\psi$ MC', r'$B \rightarrow K^{(*)} J/\psi$ data', r'')
