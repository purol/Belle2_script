import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import math

Scale_CHG_validation_MC15rd = (0.361673/1.458959)
Scale_MIX_validation_MC15rd = (0.361673/1.458959)
Scale_UUBAR_validation_MC15rd = (0.361673/1.458959)
Scale_DDBAR_validation_MC15rd = (0.361673/1.458959)
Scale_SSBAR_validation_MC15rd = (0.361673/1.458959)
Scale_CHARM_validation_MC15rd = (0.361673/1.458959)

def calculate_efficiency(data_total, data_pass, efficiency, efficiency_err):
    for N_total, N_pass in zip(data_total, data_pass):
        if (N_total != 0):
            efficiency.append(N_pass / N_total)
            efficiency_err.append((math.sqrt(N_pass * (1 - (N_pass / N_total)))) / N_total)
        else:
            efficiency.append(0.0)
            efficiency_err.append(0.0)

def calculate_efficiency_without_err(data_total, data_pass, efficiency):
    for N_total, N_pass in zip(data_total, data_pass):
        if (N_total != 0):
            efficiency.append(N_pass / N_total)
        else:
            efficiency.append(0.0)

# Define a function to plot and save a dataset
def plot_and_save_data(x, SIGNAL_nominal, Jpsi_data, Jpsi_data_err, embedded_data, embedded_data_err, filename, title):
    fig, ax = plt.subplots(figsize=(8, 7))

    trans1 = ax.transData + ScaledTranslation(-5/72, 0, fig.dpi_scale_trans)
    trans2 = ax.transData + ScaledTranslation(+5/72, 0, fig.dpi_scale_trans)
    er1 = ax.errorbar(x, Jpsi_data, yerr=Jpsi_data_err, marker="o", linestyle="none", transform=trans1, label=r"$B \rightarrow X_{s} J\psi$ data")
    er2 = ax.errorbar(x, embedded_data, yerr=embedded_data_err, marker="o", linestyle="none", transform=trans2, label=r"embedded $B \rightarrow X_{s} J\psi$ data")

    # Add horizontal lines for SIGNAL_nominal at each bin
    for i, val in enumerate(SIGNAL_nominal):
        binmin = (2.0 * i) / (2.0 * len(SIGNAL_nominal))
        binmax = binmin + (2.0 / (2.0 * len(SIGNAL_nominal)))
        ax.axhline(y=val, color='gray', linestyle='--', linewidth=0.5, xmin=binmin, xmax=binmax)

    # Create a separate line for SIGNAL_nominal values without marker
    signal_nominal_line = ax.plot([], [], color='gray', linestyle='--', linewidth=0.5, label=r'$B \rightarrow X_{s} \nu \bar{\nu}$ MC')
    
    # Rotate the x-axis labels by 45 degrees
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=30, ha='right')

    # Add a title and labels for the axes
    ax.set_xlabel('Cut criteria')
    ax.set_ylabel('efficiency')

    ax.legend(loc='best')

    ax.text(0.02, 0.95, title, transform=ax.transAxes, color='black', fontsize=16)

    plt.tight_layout()

    # Save the figure as a PNG image
    plt.savefig(filename, dpi=300)
    #plt.show()
    plt.close()  # Close the figure to release resources

x = [r'$M_{bc}^{tag} > 5.27$ GeV', r'$|\Delta E^{tag}| < 0.2$ GeV', r'$E_{ecl} < 1.3$ GeV', r'$0.297 < \theta_{miss} < 2.618$', r'$0.5 < p_{sig} < 2.96$ GeV', r'D veto']

Jpsi_data_initial = [175, 185, 75, 17, 2, 0]
Jpsi_data_Mbc = [152, 149, 47, 9, 2, 0]
Jpsi_data_deltaE = [133, 125, 39, 9, 0, 0]
Jpsi_data_Eecl = [133, 125, 39, 9, 0, 0]
Jpsi_data_missingtheta = [126, 122, 35, 9, 0, 0]
Jpsi_data_psig = [126, 122, 35, 9, 0, 0]
Jpsi_data_Dveto = [126, 122, 34, 8, 0, 0]

Jpsi_data_Mbc_eff = []
Jpsi_data_deltaE_eff = []
Jpsi_data_Eecl_eff = []
Jpsi_data_missingtheta_eff = []
Jpsi_data_psig_eff = []
Jpsi_data_Dveto_eff = []

Jpsi_data_Mbc_eff_err = []
Jpsi_data_deltaE_eff_err = []
Jpsi_data_Eecl_eff_err = []
Jpsi_data_missingtheta_eff_err = []
Jpsi_data_psig_eff_err = []
Jpsi_data_Dveto_eff_err = []

Jpsi_data_Kpi_eff = []
Jpsi_data_Kpi_eff_err = []

calculate_efficiency(Jpsi_data_initial, Jpsi_data_Mbc, Jpsi_data_Mbc_eff, Jpsi_data_Mbc_eff_err)
calculate_efficiency(Jpsi_data_Mbc, Jpsi_data_deltaE, Jpsi_data_deltaE_eff, Jpsi_data_deltaE_eff_err)
calculate_efficiency(Jpsi_data_deltaE, Jpsi_data_Eecl, Jpsi_data_Eecl_eff, Jpsi_data_Eecl_eff_err)
calculate_efficiency(Jpsi_data_Eecl, Jpsi_data_missingtheta, Jpsi_data_missingtheta_eff, Jpsi_data_missingtheta_eff_err)
calculate_efficiency(Jpsi_data_missingtheta, Jpsi_data_psig, Jpsi_data_psig_eff, Jpsi_data_psig_eff_err)
calculate_efficiency(Jpsi_data_psig, Jpsi_data_Dveto, Jpsi_data_Dveto_eff, Jpsi_data_Dveto_eff_err)

Jpsi_data_Kpi_eff.append(Jpsi_data_Mbc_eff[1])
Jpsi_data_Kpi_eff.append(Jpsi_data_deltaE_eff[1])
Jpsi_data_Kpi_eff.append(Jpsi_data_Eecl_eff[1])
Jpsi_data_Kpi_eff.append(Jpsi_data_missingtheta_eff[1])
Jpsi_data_Kpi_eff.append(Jpsi_data_psig_eff[1])
Jpsi_data_Kpi_eff.append(Jpsi_data_Dveto_eff[1])

Jpsi_data_Kpi_eff_err.append(Jpsi_data_Mbc_eff_err[1])
Jpsi_data_Kpi_eff_err.append(Jpsi_data_deltaE_eff_err[1])
Jpsi_data_Kpi_eff_err.append(Jpsi_data_Eecl_eff_err[1])
Jpsi_data_Kpi_eff_err.append(Jpsi_data_missingtheta_eff_err[1])
Jpsi_data_Kpi_eff_err.append(Jpsi_data_psig_eff_err[1])
Jpsi_data_Kpi_eff_err.append(Jpsi_data_Dveto_eff_err[1])

embedded_data_initial = [170, 185, 54, 13, 0, 1]
embedded_data_Mbc = [150, 144, 33, 5, 0, 1]
embedded_data_deltaE = [133, 121, 28, 3, 0, 1]
embedded_data_Eecl = [133, 121, 28, 3, 0, 1]
embedded_data_missingtheta = [101, 81, 16, 1, 0, 1]
embedded_data_psig = [98, 75, 15, 1, 0, 1]
embedded_data_Dveto = [98, 73, 15, 1, 0, 1]

embedded_data_Mbc_eff = []
embedded_data_deltaE_eff = []
embedded_data_Eecl_eff = []
embedded_data_missingtheta_eff = []
embedded_data_psig_eff = []
embedded_data_Dveto_eff = []

embedded_data_Mbc_eff_err = []
embedded_data_deltaE_eff_err = []
embedded_data_Eecl_eff_err = []
embedded_data_missingtheta_eff_err = []
embedded_data_psig_eff_err = []
embedded_data_Dveto_eff_err = []

embedded_data_Kpi_eff = []
embedded_data_Kpi_eff_err = []

calculate_efficiency(embedded_data_initial, embedded_data_Mbc, embedded_data_Mbc_eff, embedded_data_Mbc_eff_err)
calculate_efficiency(embedded_data_Mbc, embedded_data_deltaE, embedded_data_deltaE_eff, embedded_data_deltaE_eff_err)
calculate_efficiency(embedded_data_deltaE, embedded_data_Eecl, embedded_data_Eecl_eff, embedded_data_Eecl_eff_err)
calculate_efficiency(embedded_data_Eecl, embedded_data_missingtheta, embedded_data_missingtheta_eff, embedded_data_missingtheta_eff_err)
calculate_efficiency(embedded_data_missingtheta, embedded_data_psig, embedded_data_psig_eff, embedded_data_psig_eff_err)
calculate_efficiency(embedded_data_psig, embedded_data_Dveto, embedded_data_Dveto_eff, embedded_data_Dveto_eff_err)

embedded_data_Kpi_eff.append(embedded_data_Mbc_eff[1])
embedded_data_Kpi_eff.append(embedded_data_deltaE_eff[1])
embedded_data_Kpi_eff.append(embedded_data_Eecl_eff[1])
embedded_data_Kpi_eff.append(embedded_data_missingtheta_eff[1])
embedded_data_Kpi_eff.append(embedded_data_psig_eff[1])
embedded_data_Kpi_eff.append(embedded_data_Dveto_eff[1])

embedded_data_Kpi_eff_err.append(embedded_data_Mbc_eff_err[1])
embedded_data_Kpi_eff_err.append(embedded_data_deltaE_eff_err[1])
embedded_data_Kpi_eff_err.append(embedded_data_Eecl_eff_err[1])
embedded_data_Kpi_eff_err.append(embedded_data_missingtheta_eff_err[1])
embedded_data_Kpi_eff_err.append(embedded_data_psig_eff_err[1])
embedded_data_Kpi_eff_err.append(embedded_data_Dveto_eff_err[1])

Signal_MC_initial = [69.56361, 60.85152, 26.31493, 4.525909, 0.364424, 0.347064]
Signal_MC_Mbc = [52.05795, 40.56013, 14.9245, 2.290768, 0.190013, 0.247485]
Signal_MC_deltaE = [47.28268, 36.43381, 13.26239, 2.002802, 0.171587, 0.224066]
Signal_MC_Eecl = [45.43815, 35.50649, 12.98649, 1.96625, 0.170176, 0.223127]
Signal_MC_missingtheta = [41.89615, 33.62171, 12.61645, 1.933378, 0.169394, 0.218523]
Signal_MC_psig = [39.95194, 32.85438, 12.4467, 1.895295, 0.168522, 0.216476]
Signal_MC_Dveto = [39.95194, 32.76187, 12.13724, 1.756578, 0.15261, 0.199373]

Signal_MC_Mbc_eff = []
Signal_MC_deltaE_eff = []
Signal_MC_Eecl_eff = []
Signal_MC_missingtheta_eff = []
Signal_MC_psig_eff = []
Signal_MC_Dveto_eff = []

Signal_MC_Kpi_eff = []

calculate_efficiency_without_err(Signal_MC_initial, Signal_MC_Mbc, Signal_MC_Mbc_eff)
calculate_efficiency_without_err(Signal_MC_Mbc, Signal_MC_deltaE, Signal_MC_deltaE_eff)
calculate_efficiency_without_err(Signal_MC_deltaE, Signal_MC_Eecl, Signal_MC_Eecl_eff)
calculate_efficiency_without_err(Signal_MC_Eecl, Signal_MC_missingtheta, Signal_MC_missingtheta_eff)
calculate_efficiency_without_err(Signal_MC_missingtheta, Signal_MC_psig, Signal_MC_psig_eff)
calculate_efficiency_without_err(Signal_MC_psig, Signal_MC_Dveto, Signal_MC_Dveto_eff)

Signal_MC_Kpi_eff.append(Signal_MC_Mbc_eff[1])
Signal_MC_Kpi_eff.append(Signal_MC_deltaE_eff[1])
Signal_MC_Kpi_eff.append(Signal_MC_Eecl_eff[1])
Signal_MC_Kpi_eff.append(Signal_MC_missingtheta_eff[1])
Signal_MC_Kpi_eff.append(Signal_MC_psig_eff[1])
Signal_MC_Kpi_eff.append(Signal_MC_Dveto_eff[1])

plot_and_save_data(x, Signal_MC_Kpi_eff, Jpsi_data_Kpi_eff, Jpsi_data_Kpi_eff_err, embedded_data_Kpi_eff, embedded_data_Kpi_eff_err, "eff_flow_compare.png", "")


#plot_and_save_data(x, y1_1, yerr1_1, "FBDT_eff.png", r'$B \rightarrow X_{s} J/\psi$', r'')
