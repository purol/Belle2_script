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
    ax.set_ylim(-1, y_max * 1.2)
    ax.set_xlim(-0.5, len(x) - 0.5)
    
    # Rotate the x-axis labels by 45 degrees
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=30, ha='right', fontsize=14)

    # Add a title and labels for the axes
    ax.set_xlabel('Decay Mode', fontsize=14)
    ax.set_ylabel('Fraction at 7 decays[%]', fontsize=14)

    ax.legend(loc='best', fontsize=15)

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

x = [r'$K 1\pi$ w/o $\pi^{0}$', r'$K 1\pi$ w/ $1\pi^{0}$', r'$K 2\pi$ w/o $\pi^{0}$', r'$K 2\pi$ w/ $1\pi^{0}$', r'$K 3\pi$ w/o $\pi^{0}$', r'$K 3\pi$ w/ $1\pi^{0}$', r'$K 4\pi$ at most $1\pi^{0}$']

y1_1 = [11.2, 6.18, 24.4, 46.7, 0.49, 10.5, 0.55] # bn1480_v4 [1.15, 1.5]
yerr1_1 = [0.73, 0.61, 1.35, 1.87, 0.31, 1.23, 0.39]
y1_2 = [3.79, 1.4, 19.1, 25.1, 11.5, 32.7, 6.43] # bn1480_v4 [1.5, 2.0]
yerr1_2 = [0.40, 0.4, 1.31, 1.86, 1.09, 2.40, 1.78]


y2_1 = [12.5, 7.4, 29.5, 29.7, 4.7, 14.5, 1.5] # Xs Jpsi MC + correction [1.15, 1.5]
yerr2_1 = [0.5, 0.7, 1.0, 2.4, 0.6, 2.6, 0.0]
y2_2 = [5.5, 3.3, 17.9, 19.0, 12.6, 27.0, 12.5] # Xs Jpsi MC + correction [1.5, 2.0]
yerr2_2 = [0.20, 0.31, 0.61, 1.5, 1.5, 4.9, 0.0]


plot_and_save_data(x, y1_1, yerr1_1, y2_1, yerr2_1, "fraction_1.15_1.5.png", r'$B \rightarrow X_{s} \gamma$ [bn1480_v4]', r'$B \rightarrow X_{s} J/\psi$ MC + correction', r'$1.15 < M_{X_{s}} < 1.5$ GeV')
plot_and_save_data(x, y1_2, yerr1_2, y2_2, yerr2_2, "fraction_1.5_2.0.png", r'$B \rightarrow X_{s} \gamma$ [bn1480_v4]', r'$B \rightarrow X_{s} J/\psi$ + correction', r'$1.5 < M_{X_{s}} < 2.0$ GeV')

print_ratio(y1_1, yerr1_1, y2_1, "weight_fragmentation_1.15_1.5.txt")
print_ratio(y1_2, yerr1_2, y2_2, "weight_fragmentation_1.5_2.0.txt")
