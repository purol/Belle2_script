import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import math

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
    ax.set_ylabel('Fraction [%]')

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

x = [r'$K^{+*}$', r'$K^{+}(1270)$', r'$K^{+}\pi^{-}\pi^{+}$', r'$\eta K^{+}$', r'$\phi K^{+}$', r'$\eta^{\prime} K^{+}$']

# B->Xs gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S049.225&home=MXXX045
BR_Xsgamma = 0.000349
BR_sigma_Xsgamma = 0.000019
BR_rel_Xsgamma = BR_sigma_Xsgamma / BR_Xsgamma

# B->K*+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.9&home=MXXX045
BR_Kstargamma = 0.0000392
BR_sigma_Kstargamma = 0.0000022
BR_rel_Kstargamma = BR_sigma_Kstargamma / BR_Kstargamma

# B->K1(1270)+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.53&home=MXXX045
BR_K1270gamma = 0.000044
BR_sigma_K1270gamma = 0.0000065
BR_rel_K1270gamma = BR_sigma_K1270gamma / BR_K1270gamma

# B->K+ pi- pi+ gamma (all resonance), https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.307&home=MXXX045
BR_Kpipigamma = 0.0000258
BR_sigma_Kpipigamma = 0.0000015
BR_rel_Kpipigamma = BR_sigma_Kpipigamma / BR_Kpipigamma

# B->eta K+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.403&home=MXXX045
BR_etaKgamma = 0.0000079
BR_sigma_etaKgamma = 0.0000009
BR_rel_etaKgamma = BR_sigma_etaKgamma / BR_etaKgamma

# B->phi K+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.347&home=MXXX045
BR_phiKgamma = 0.0000027
BR_sigma_phiKgamma = 0.0000004
BR_rel_phiKgamma = BR_sigma_phiKgamma / BR_phiKgamma

# B->eta' K+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.451&home=MXXX045
BR_etaprimeKgamma = 0.0000029
BR_sigma_etaprimeKgamma = (0.000001 + 0.0000009)/2.0
BR_rel_etaprimeKgamma = BR_sigma_etaprimeKgamma / BR_etaprimeKgamma



# B=direct=>X J/psi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S049.23&home=MXXX045
BR_XJpsi = 0.0078
BR_sigma_XJpsi = 0.0004
BR_rel_XJpsi = BR_sigma_XJpsi / BR_XJpsi

# B=direct=>K*+ J/psi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.142&home=MXXX045
BR_KstarJpsi = 0.00143
BR_sigma_KstarJpsi = 0.00008
BR_rel_KstarJpsi = BR_sigma_KstarJpsi / BR_KstarJpsi

# B->K1(1270)+ Jpsi, https://pdglive.lbl.gov/Particle.action?init=0&node=S041&home=MXXX045
BR_K1270Jpsi = 0.0018
BR_sigma_K1270Jpsi = 0.0005
BR_rel_K1270Jpsi = BR_sigma_K1270Jpsi / BR_K1270Jpsi

# B=direct=>K+ pi- pi+ J/psi (all resonance), https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.19&home=MXXX045
BR_KpipiJpsi = 0.00081
BR_sigma_KpipiJpsi = 0.00013
BR_rel_KpipiJpsi = BR_sigma_KpipiJpsi / BR_KpipiJpsi

# B=direct=>eta K+ Jpsi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.354&home=MXXX045
BR_etaKJpsi = 0.000124
BR_sigma_etaKJpsi = 0.000014
BR_rel_etaKJpsi = BR_sigma_etaKJpsi / BR_etaKJpsi

# B->phi K+ Jpsi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.268&home=MXXX045
BR_phiKJpsi = 0.00005
BR_sigma_phiKJpsi = 0.000004
BR_rel_phiKJpsi = BR_sigma_phiKJpsi / BR_phiKJpsi

# B->eta' K+ Jpsi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.477&home=MXXX045
BR_etaprimeJpsi = 0.000031
BR_sigma_etaprimeKJpsi = 0.000004
BR_rel_etaprimeKJpsi = BR_sigma_etaprimeKJpsi / BR_etaprimeJpsi


y1 = [BR_Kstargamma, BR_K1270gamma, BR_Kpipigamma, BR_etaKgamma, BR_phiKgamma, BR_etaprimeKgamma] # for B->Xs gamma
y1 = [x / BR_Xsgamma for x in y1]
yrelerr1 = [math.sqrt(BR_rel_Xsgamma**2 + BR_rel_Kstargamma**2), math.sqrt(BR_rel_Xsgamma**2 + BR_rel_K1270gamma**2), math.sqrt(BR_rel_Xsgamma**2 + BR_rel_Kpipigamma**2),\
         math.sqrt(BR_rel_Xsgamma**2 + BR_rel_etaKgamma**2), math.sqrt(BR_rel_Xsgamma**2 + BR_rel_phiKgamma**2), math.sqrt(BR_rel_Xsgamma**2 + BR_rel_etaprimeKgamma**2)]
yerr1 = [x * y for x, y in zip(y1, yrelerr1)]

y2 = [BR_KstarJpsi, BR_K1270Jpsi, BR_KpipiJpsi, BR_etaKJpsi, BR_phiKJpsi, BR_etaprimeJpsi] # for B->Xs J/psi
y2 = [x / BR_XJpsi for x in y2]
yrelerr2 = [math.sqrt(BR_rel_XJpsi**2 + BR_rel_KstarJpsi**2), math.sqrt(BR_rel_XJpsi**2 + BR_rel_K1270Jpsi**2), math.sqrt(BR_rel_XJpsi**2 + BR_rel_KpipiJpsi**2),\
         math.sqrt(BR_rel_XJpsi**2 + BR_rel_etaKJpsi**2), math.sqrt(BR_rel_XJpsi**2 + BR_rel_phiKJpsi**2), math.sqrt(BR_rel_XJpsi**2 + BR_rel_etaprimeKJpsi**2)]
yerr2 = [x * y for x, y in zip(y2, yrelerr2)]

plot_and_save_data(x, y1, yerr1, y2, yerr2, "fragmentation_compare", r'$B \rightarrow X_{s} \gamma$', r'$B \rightarrow X_{s} J/\psi$', r'')
