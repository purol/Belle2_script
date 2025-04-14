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
    ax.set_xticklabels(x, rotation=30, ha='right')

    # Add a title and labels for the axes
    ax.set_xlabel('Decay Mode', fontsize=14)
    ax.set_ylabel('Fraction at 9 decays[%]', fontsize=14)

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

x = [r'$K 1\pi$ w/o $\pi^{0}$', r'$K 1\pi$ w/ $1\pi^{0}$', r'$K 2\pi$ w/o $\pi^{0}$', r'$K 2\pi$ w/ $1\pi^{0}$', r'$K 3\pi$ w/o $\pi^{0}$', r'$K 3\pi$ w/ $1\pi^{0}$', r'$K 4\pi$ at most $1\pi^{0}$', r'including $2\pi^{0}$', r'including $3K$']

y1 = [4.43, 2.31, 14.7, 21.9, 6.12, 18.6, 8.64, 19.7, 2.05] # bn1480_v4 [1.15, 2.8]
yerr1 = [0.33, 0.27, 0.95, 1.45, 0.71, 1.65, 2.07, 3.61, 0.3]
y1_1 = [10.6, 5.86, 23.1, 44.3, 0.46, 9.96, 0.52, 4.76, 0.0] # bn1480_v4 [1.15, 1.5]
yerr1_1 = [0.69, 0.58, 1.28, 1.77, 0.29, 1.17, 0.37, 1.14, 0.0]
y1_2 = [3.04, 1.1, 15.3, 20.1, 9.24, 26.2, 5.15, 15.1, 2.14] # bn1480_v4 [1.5, 2.0]
yerr1_2 = [0.32, 0.3, 1.05, 1.49, 0.87, 1.92, 1.43, 2.59, 0.25]
y1_3 = [1.6, 0.49, 8.73, 7.72, 16.3, 25.8, 9.93, 21.1, 5.21] # bn1480_v4 [2.0, 2.4]
yerr1_3 = [0.68, 0.71, 2.31, 3.09, 4.09, 6.89, 8.35, 14.1, 1.51]

y2 = [39.00, 19.48, 16.61, 17.04, 0.40, 3.03, 0.05, 4.40, 0.00] # sigal MC [1.1, inf]
yerr2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
y2_1 = [30.97, 15.49, 17.38, 20.90, 1.38, 8.15, 0.75, 4.90, 0.02] # sigal MC [1.15, 1.5]
yerr2_1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
y2_2 = [20.40, 10.20, 15.18, 19.65, 3.63, 13.90, 5.20, 8.54, 3.31] # sigal MC [1.5, 2.0]
yerr2_2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
y2_3 = [12.75, 6.37, 12.14, 15.17, 5.89, 15.95, 12.49, 14.19, 5.07] # sigal MC [2.0, 2.4]
yerr2_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

plot_and_save_data(x, y1_1, yerr1_1, y2_1, yerr2_1, "fraction_1.15_1.5.png", r'$B \rightarrow X_{s} \gamma$ [bn1480_v4]', r'$B \rightarrow X_{s} \nu \bar{\nu}$ signal MC', r'$1.15 < M_{X_{s}} < 1.5$ GeV')
plot_and_save_data(x, y1_2, yerr1_2, y2_2, yerr2_2, "fraction_1.5_2.0.png", r'$B \rightarrow X_{s} \gamma$ [bn1480_v4]', r'$B \rightarrow X_{s} \nu \bar{\nu}$ signal MC', r'$1.5 < M_{X_{s}} < 2.0$ GeV')
plot_and_save_data(x, y1_3, yerr1_3, y2_3, yerr2_3, "fraction_2.0_2.4.png", r'$B \rightarrow X_{s} \gamma$ [bn1480_v4]', r'$B \rightarrow X_{s} \nu \bar{\nu}$ signal MC', r'$2.0 < M_{X_{s}} < 2.4$ GeV')

print_ratio(y1_1, yerr1_1, y2_1, "weight_fragmentation_1.15_1.5.txt")
print_ratio(y1_2, yerr1_2, y2_2, "weight_fragmentation_1.5_2.0.txt")
print_ratio(y1_3, yerr1_3, y2_3, "weight_fragmentation_2.0_2.4.txt")
