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

x = [r'$K^{+*}$', r'$K^{+}(1270)$', r'$K^{+}\pi^{-}\pi^{+}$', r'$\eta K^{+}$']

# B->Xs gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S049.225&home=MXXX045
BR_Xsgamma = 0.000349
BR_sigma_Xsgamma = 0.000019

# B->K*+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.9&home=MXXX045
BR_Kstargamma = 0.0000392
BR_sigma_Kstargamma = 0.0000022

# B->K1(1270)+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.53&home=MXXX045
BR_K1270gamma = 0.000044
BR_sigma_K1270gamma = 0.0000065

# B->K+ pi- pi+ gamma (all resonance), https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.307&home=MXXX045
BR_Kpipigamma = 0.0000258
BR_sigma_Kpipigamma = 0.0000015

# B->eta K+ gamma, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.403&home=MXXX045
BR_etaKgamma = 0.0000079
BR_sigma_etaKgamma = 0.0000009



# B=direct=>X J/psi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S049.23&home=MXXX045
BR_XJpsi = 0.0078
BR_sigma_XJpsi = 0.0004

# B=direct=>K*+ J/psi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.142&home=MXXX045
BR_KstarJpsi = 0.00143
BR_sigma_KstarJpsi = 0.00008

# B->K1(1270)+ Jpsi, https://pdglive.lbl.gov/Particle.action?init=0&node=S041&home=MXXX045
BR_K1270Jpsi = 0.0018
BR_sigma_K1270Jpsi =0.0005

# B=direct=>K+ pi- pi+ J/psi (all resonance), https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.19&home=MXXX045
BR_KpipiJpsi = 0.00081
BR_sigma_KpipiJpsi = 0.00013

# B=direct=>eta K+ Jpsi, https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S041.354&home=MXXX045
BR_etaKJpsi = 0.000124
BR_sigma_etaKJpsi = 0.000014


# BR(K1+(1270) -> K*0 pi+) * BR(K*0 -> K+ pi-)
BR_K1270toKstrapi_BRK0startoKpluspiminus = (0.21 * (2.0/3.0)) * (0.99754 * (2.0/3.0))
BR_sigma_K1270toKstrapi_BRK0startoKpluspiminus = math.sqrt((10.0/21.0)**2 + (0.021/99.754)**2) * BR_K1270toKstrapi_BRK0startoKpluspiminus

# BR(K1+(1270) -> K+ rho0) * BR(rho0 -> pi+ pi-)
BR_K1270toKplusrhozero_BRrhozerotopipluspiminus = (0.38 * (1.0/3.0)) * (1.0)
BR_sigma_K1270toKplusrhozero_BRrhozerotopipluspiminus = (13.0/38.0) * BR_K1270toKplusrhozero_BRrhozerotopipluspiminus

# BR(K1+(1270) -> K+ omega0) * BR(omega0 -> pi+ pi-)
BR_K1270toKplusomegazero_BRomegazerotopipluspiminus = (0.11) * (0.0153)
BR_sigma_K1270toKplusomegazero_BRomegazerotopipluspiminus = math.sqrt((2.0/11.0)**2 + (0.12/1.53)**2) * BR_K1270toKplusomegazero_BRomegazerotopipluspiminus

# BR(K1+(1270) -> K0*(1430) pi+) * BR(K0*(1430) -> K+ pi-)
BR_K1270toKzerostar1430piplus_BRKzerostar1430toKpluspiminus = (0.28 * (2.0/3.0)) * (0.93 * (2.0/3.0))
BR_sigma_K1270toKzerostar1430piplus_BRKzerostar1430toKpluspiminus = math.sqrt((0.04/0.28)**2 + (0.1/0.93)**2) * BR_K1270toKzerostar1430piplus_BRKzerostar1430toKpluspiminus

Total_Decayfraction = 1.40736
DecayFraction_1 = 0.232
DecayFraction_1_sigma = math.sqrt(0.017**2+0.058**2)
DecayFraction_2 = 0.383
DecayFraction_2_sigma = math.sqrt(0.016**2+0.036**2)
DecayFraction_3 = 0.0045
DecayFraction_3_sigma = math.sqrt(0.0017**2+0.0014**2)
DecayFraction_4 = 0.0157
DecayFraction_4_sigma = math.sqrt(0.0052**2+0.0049**2)

BR_K1270Jpsi_1 = BR_KpipiJpsi * (DecayFraction_1 / Total_Decayfraction) / BR_K1270toKstrapi_BRK0startoKpluspiminus
BR_K1270Jpsi_1_sigma = math.sqrt((BR_sigma_KpipiJpsi/BR_KpipiJpsi)**2 + (DecayFraction_1_sigma/DecayFraction_1)**2 + (BR_sigma_K1270toKstrapi_BRK0startoKpluspiminus/BR_K1270toKstrapi_BRK0startoKpluspiminus)**2 ) * BR_K1270Jpsi_1

BR_K1270Jpsi_2 = BR_KpipiJpsi * (DecayFraction_2 / Total_Decayfraction) / BR_K1270toKplusrhozero_BRrhozerotopipluspiminus
BR_K1270Jpsi_2_sigma = math.sqrt((BR_sigma_KpipiJpsi/BR_KpipiJpsi)**2 + (DecayFraction_2_sigma/DecayFraction_2)**2 + (BR_sigma_K1270toKplusrhozero_BRrhozerotopipluspiminus/BR_K1270toKplusrhozero_BRrhozerotopipluspiminus)**2 ) * BR_K1270Jpsi_2

BR_K1270Jpsi_3 = BR_KpipiJpsi * (DecayFraction_3 / Total_Decayfraction) / BR_K1270toKplusomegazero_BRomegazerotopipluspiminus
BR_K1270Jpsi_3_sigma = math.sqrt((BR_sigma_KpipiJpsi/BR_KpipiJpsi)**2 + (DecayFraction_3_sigma/DecayFraction_3)**2 + (BR_sigma_K1270toKplusomegazero_BRomegazerotopipluspiminus/BR_K1270toKplusomegazero_BRomegazerotopipluspiminus)**2 ) * BR_K1270Jpsi_3

BR_K1270Jpsi_4 = BR_KpipiJpsi * (DecayFraction_4 / Total_Decayfraction) / BR_K1270toKzerostar1430piplus_BRKzerostar1430toKpluspiminus
BR_K1270Jpsi_4_sigma = math.sqrt((BR_sigma_KpipiJpsi/BR_KpipiJpsi)**2 + (DecayFraction_4_sigma/DecayFraction_4)**2 + (BR_sigma_K1270toKzerostar1430piplus_BRKzerostar1430toKpluspiminus/BR_K1270toKzerostar1430piplus_BRKzerostar1430toKpluspiminus)**2 ) * BR_K1270Jpsi_4

print(BR_K1270Jpsi_1)
print(BR_K1270Jpsi_1_sigma)

print(BR_K1270Jpsi_2)
print(BR_K1270Jpsi_2_sigma)

print(BR_K1270Jpsi_3)
print(BR_K1270Jpsi_3_sigma)

print(BR_K1270Jpsi_4)
print(BR_K1270Jpsi_4_sigma)


