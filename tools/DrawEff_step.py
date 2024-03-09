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

def relative_err_if_data(num):
    if (num == 0):
        return 0
    elif(num > 0):
        return (math.sqrt(num) / num)
    else:
        print("Nevt is negative!")
        exit(1)

def efficiency(Nevt_after, Nevt_before):
    if(Nevt_before == 0):
        return 0
    else:
        return (Nevt_after/Nevt_before)

# Define a function to plot and save a dataset
def plot_and_save_data(x, SIGNAL_nominal, Jpsi_data, Jpsi_data_err, Jpsi_MC, Jpsi_MC_err, embedded_data, embedded_data_err, embedded_MC, embedded_MC_Err, filename, title):
    fig, ax = plt.subplots(figsize=(8, 7))

    trans1 = ax.transData + ScaledTranslation(-10/72, 0, fig.dpi_scale_trans)
    trans2 = ax.transData + ScaledTranslation(-5/72, 0, fig.dpi_scale_trans)
    trans3 = ax.transData + ScaledTranslation(+5/72, 0, fig.dpi_scale_trans)
    trans4 = ax.transData + ScaledTranslation(+10/72, 0, fig.dpi_scale_trans)
    er1 = ax.errorbar(x, Jpsi_data, yerr=Jpsi_data_err, marker="o", linestyle="none", transform=trans1, label=r"$B \rightarrow X_{s} J\psi$ data")
    er2 = ax.errorbar(x, Jpsi_MC, yerr=Jpsi_MC_err, marker="o", linestyle="none", transform=trans2, label=r"$B \rightarrow X_{s} J\psi$ MC")
    er3 = ax.errorbar(x, embedded_data, yerr=embedded_data_err, marker="o", linestyle="none", transform=trans3, label=r"embedded $B \rightarrow X_{s} J\psi$ data")
    er4 = ax.errorbar(x, embedded_MC, yerr=embedded_MC_Err, marker="o", linestyle="none", transform=trans4, label=r"embedded $B \rightarrow X_{s} J\psi$ MC")
    
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

Jpsi_data_initial_eff = []
Jpsi_data_Mbc_eff = []
Jpsi_data_deltaE_eff = []
Jpsi_data_Eecl_eff = []
Jpsi_data_missingtheta_eff = []
Jpsi_data_psig_eff = []
Jpsi_data_Dveto_eff = []

embedded_data_initial = [170, 185, 54, 13, 0, 1]
embedded_data_Mbc = [150, 144, 33, 5, 0, 1]
embedded_data_deltaE = [133, 121, 28, 3, 0, 1]
embedded_data_Eecl = [133, 121, 28, 3, 0, 1]
embedded_data_missingtheta = [101, 81, 16, 1, 0, 1]
embedded_data_psig = [98, 75, 15, 1, 0, 1]
embedded_data_Dveto = [98, 73, 15, 1, 0, 1]

embedded_data_initial_eff = []
embedded_data_Mbc_eff = []
embedded_data_deltaE_eff = []
embedded_data_Eecl_eff = []
embedded_data_missingtheta_eff = []
embedded_data_psig_eff = []
embedded_data_Dveto_eff = []

# nominal MC
SIGNAL_nominal_MC_before = [24.50669, 19.15963, 6.137504, 0.801534, 0.046628, 0.118388]
SIGNAL_nominal_MC_after = [13.74668, 6.419036, 1.890514, 0.219984, 0.009182, 0.040427]

# J/psi data
Jpsi_data_before = [79, 75, 20, 6, 0, 0]
Jpsi_data_after = [58, 30, 4, 1, 0, 0]

# J/psi MC
Jpsi_SIGNAL_MC_before = [125.6843, 99.90289, 59.49552, 3.71847, 0, 0.743694]
Jpsi_CHG_MC_before = [0.247898, 0.495796, 0, 0, 0, 0]
Jpsi_MIX_MC_before = [0, 0, 0.743694, 0, 0, 0]
Jpsi_UUBAR_MC_before = [0.247898, 0.247898, 0.991592, 0.247898, 0, 0]
Jpsi_DDBAR_MC_before = [0, 0, 0.247898, 0, 0, 0]
Jpsi_SSBAR_MC_before = [0, 0.247898, 0.247898, 0, 0, 0]
Jpsi_CHARM_MC_before = [0, 0, 0.743694, 0, 0, 0]
Jpsi_SIGNAL_MC_after = [97.67181, 32.97043, 13.13859, 0.495796, 0, 0]
Jpsi_CHG_MC_after = [0, 0.247898, 0, 0, 0, 0]
Jpsi_MIX_MC_after = [0, 0, 0, 0, 0, 0]
Jpsi_UUBAR_MC_after = [0, 0, 0, 0, 0, 0]
Jpsi_DDBAR_MC_after = [0, 0, 0, 0, 0, 0]
Jpsi_SSBAR_MC_after = [0, 0, 0, 0, 0, 0]
Jpsi_CHARM_MC_after = [0, 0, 0, 0, 0, 0]

# embedded data
embedded_data_before = [56, 47, 7, 1, 0, 1]
embedded_data_after = [36, 5, 0, 0, 0, 0]

# embedded MC
embedded_CHG_MC_before = [86.2685, 41.64686, 9.668022, 0.743694, 0, 0.495796]
embedded_MIX_MC_before = [5.205858, 32.72253, 4.710062, 0.495796, 0, 0]
embedded_UUBAR_MC_before = [0.743694, 0.991592, 0, 0, 0, 0]
embedded_DDBAR_MC_before = [0.247898, 0.247898, 0, 0, 0, 0]
embedded_SSBAR_MC_before = [0.247898, 0.247898, 0, 0, 0, 0]
embedded_CHARM_MC_before = [0.247898, 0.247898, 0.247898, 0, 0, 0]
embedded_CHG_MC_after = [55.03335, 3.966368, 0.247898, 0, 0, 0]
embedded_MIX_MC_after = [1.735286, 1.23949, 0, 0, 0, 0]
embedded_UUBAR_MC_after = [0.495796, 0.247898, 0, 0, 0, 0]
embedded_DDBAR_MC_after = [0, 0, 0, 0, 0, 0]
embedded_SSBAR_MC_after = [0, 0, 0, 0, 0, 0]
embedded_CHARM_MC_after = [0, 0, 0, 0, 0, 0]


# nominal MC
SIGNAL_nominal_MC_eff = []

# J/psi data
Jpsi_data_eff = []
Jpsi_data_eff_err = []

# J/psi MC
Jpsi_MC_eff = []
Jpsi_MC_eff_err = []

# embedded data
embedded_data_eff = []
embedded_data_eff_err = []

# embedded MC
embedded_MC_eff = []
embedded_MC_eff_err = []


# calculate efficiency and its uncertainty
for N_SIGNAL_nominal_MC_before, N_SIGNAL_nominal_MC_after in zip(SIGNAL_nominal_MC_before, SIGNAL_nominal_MC_after):
    SIGNAL_nominal_MC_eff.append(efficiency(N_SIGNAL_nominal_MC_after, N_SIGNAL_nominal_MC_before))

for N_Jpsi_data_before, N_Jpsi_data_after in zip(Jpsi_data_before, Jpsi_data_after):
    N_Jpsi_data_before_relative_err = relative_err_if_data(N_Jpsi_data_before)
    N_Jpsi_data_after_relative_err = relative_err_if_data(N_Jpsi_data_after)
    Jpsi_data_eff.append(efficiency(N_Jpsi_data_after, N_Jpsi_data_before))
    Jpsi_data_eff_err.append( \
        math.sqrt(N_Jpsi_data_before_relative_err*N_Jpsi_data_before_relative_err+N_Jpsi_data_after_relative_err*N_Jpsi_data_after_relative_err) * \
        efficiency(N_Jpsi_data_after, N_Jpsi_data_before)
        )
        
for N_Jpsi_SIGNAL_MC_before, \
    N_Jpsi_CHG_MC_before, \
    N_Jpsi_MIX_MC_before, \
    N_Jpsi_UUBAR_MC_before, \
    N_Jpsi_DDBAR_MC_before, \
    N_Jpsi_SSBAR_MC_before, \
    N_Jpsi_CHARM_MC_before, \
    N_Jpsi_SIGNAL_MC_after, \
    N_Jpsi_CHG_MC_after, \
    N_Jpsi_MIX_MC_after, \
    N_Jpsi_UUBAR_MC_after, \
    N_Jpsi_DDBAR_MC_after, \
    N_Jpsi_SSBAR_MC_after, \
    N_Jpsi_CHARM_MC_after \
    in zip(Jpsi_SIGNAL_MC_before,
           Jpsi_CHG_MC_before, \
           Jpsi_MIX_MC_before, \
           Jpsi_UUBAR_MC_before, \
           Jpsi_DDBAR_MC_before, \
           Jpsi_SSBAR_MC_before, \
           Jpsi_CHARM_MC_before, \
           Jpsi_SIGNAL_MC_after, \
           Jpsi_CHG_MC_after, \
           Jpsi_MIX_MC_after, \
           Jpsi_UUBAR_MC_after, \
           Jpsi_DDBAR_MC_after, \
           Jpsi_SSBAR_MC_after, \
           Jpsi_CHARM_MC_after):
    N_Jpsi_SIGNAL_MC_before_relative_err = relative_err_if_data(N_Jpsi_SIGNAL_MC_before) * (Scale_CHG_validation_MC15rd + Scale_MIX_validation_MC15rd) / 2.0
    N_Jpsi_CHG_MC_before_relative_err = relative_err_if_data(N_Jpsi_CHG_MC_before) * Scale_CHG_validation_MC15rd
    N_Jpsi_MIX_MC_before_relative_err = relative_err_if_data(N_Jpsi_MIX_MC_before) * Scale_MIX_validation_MC15rd
    N_Jpsi_UUBAR_MC_before_relative_err = relative_err_if_data(N_Jpsi_UUBAR_MC_before) * Scale_UUBAR_validation_MC15rd
    N_Jpsi_DDBAR_MC_before_relative_err = relative_err_if_data(N_Jpsi_DDBAR_MC_before) * Scale_DDBAR_validation_MC15rd
    N_Jpsi_SSBAR_MC_before_relative_err = relative_err_if_data(N_Jpsi_SSBAR_MC_before) * Scale_SSBAR_validation_MC15rd
    N_Jpsi_CHARM_MC_before_relative_err = relative_err_if_data(N_Jpsi_CHARM_MC_before) * Scale_CHARM_validation_MC15rd
    N_Jpsi_MC_before = \
                     N_Jpsi_SIGNAL_MC_before + \
                     N_Jpsi_CHG_MC_before + \
                     N_Jpsi_MIX_MC_before + \
                     N_Jpsi_UUBAR_MC_before + \
                     N_Jpsi_DDBAR_MC_before + \
                     N_Jpsi_SSBAR_MC_before + \
                     N_Jpsi_CHARM_MC_before
    if(N_Jpsi_MC_before != 0):
        N_Jpsi_MC_before_relative_err = (math.sqrt( \
            N_Jpsi_SIGNAL_MC_before_relative_err * N_Jpsi_SIGNAL_MC_before_relative_err * N_Jpsi_SIGNAL_MC_before * N_Jpsi_SIGNAL_MC_before + \
            N_Jpsi_CHG_MC_before_relative_err * N_Jpsi_CHG_MC_before_relative_err * N_Jpsi_CHG_MC_before * N_Jpsi_CHG_MC_before + \
            N_Jpsi_MIX_MC_before_relative_err * N_Jpsi_MIX_MC_before_relative_err * N_Jpsi_MIX_MC_before * N_Jpsi_MIX_MC_before + \
            N_Jpsi_UUBAR_MC_before_relative_err * N_Jpsi_UUBAR_MC_before_relative_err * N_Jpsi_UUBAR_MC_before * N_Jpsi_UUBAR_MC_before + \
            N_Jpsi_DDBAR_MC_before_relative_err * N_Jpsi_DDBAR_MC_before_relative_err * N_Jpsi_DDBAR_MC_before * N_Jpsi_DDBAR_MC_before + \
            N_Jpsi_SSBAR_MC_before_relative_err * N_Jpsi_SSBAR_MC_before_relative_err * N_Jpsi_SSBAR_MC_before * N_Jpsi_SSBAR_MC_before + \
            N_Jpsi_CHARM_MC_before_relative_err * N_Jpsi_CHARM_MC_before_relative_err * N_Jpsi_CHARM_MC_before * N_Jpsi_CHARM_MC_before \
            )) / (N_Jpsi_MC_before)
    else:
         N_Jpsi_MC_before_relative_err = 0
    N_Jpsi_SIGNAL_MC_after_relative_err = relative_err_if_data(N_Jpsi_SIGNAL_MC_after) * (Scale_CHG_validation_MC15rd + Scale_MIX_validation_MC15rd) / 2.0
    N_Jpsi_CHG_MC_after_relative_err = relative_err_if_data(N_Jpsi_CHG_MC_after) * Scale_CHG_validation_MC15rd
    N_Jpsi_MIX_MC_after_relative_err = relative_err_if_data(N_Jpsi_MIX_MC_after) * Scale_MIX_validation_MC15rd
    N_Jpsi_UUBAR_MC_after_relative_err = relative_err_if_data(N_Jpsi_UUBAR_MC_after) * Scale_UUBAR_validation_MC15rd
    N_Jpsi_DDBAR_MC_after_relative_err = relative_err_if_data(N_Jpsi_DDBAR_MC_after) * Scale_DDBAR_validation_MC15rd
    N_Jpsi_SSBAR_MC_after_relative_err = relative_err_if_data(N_Jpsi_SSBAR_MC_after) * Scale_SSBAR_validation_MC15rd
    N_Jpsi_CHARM_MC_after_relative_err = relative_err_if_data(N_Jpsi_CHARM_MC_after) * Scale_CHARM_validation_MC15rd
    N_Jpsi_MC_after = \
        N_Jpsi_SIGNAL_MC_after + \
        N_Jpsi_CHG_MC_after + \
        N_Jpsi_MIX_MC_after + \
        N_Jpsi_UUBAR_MC_after + \
        N_Jpsi_DDBAR_MC_after + \
        N_Jpsi_SSBAR_MC_after + \
        N_Jpsi_CHARM_MC_after
    if(N_Jpsi_MC_after != 0):
        N_Jpsi_MC_after_relative_err = (math.sqrt( \
            N_Jpsi_SIGNAL_MC_after_relative_err * N_Jpsi_SIGNAL_MC_after_relative_err * N_Jpsi_SIGNAL_MC_after * N_Jpsi_SIGNAL_MC_after + \
            N_Jpsi_CHG_MC_after_relative_err * N_Jpsi_CHG_MC_after_relative_err * N_Jpsi_CHG_MC_after * N_Jpsi_CHG_MC_after + \
            N_Jpsi_MIX_MC_after_relative_err * N_Jpsi_MIX_MC_after_relative_err * N_Jpsi_MIX_MC_after * N_Jpsi_MIX_MC_after + \
            N_Jpsi_UUBAR_MC_after_relative_err * N_Jpsi_UUBAR_MC_after_relative_err * N_Jpsi_UUBAR_MC_after * N_Jpsi_UUBAR_MC_after + \
            N_Jpsi_DDBAR_MC_after_relative_err * N_Jpsi_DDBAR_MC_after_relative_err * N_Jpsi_DDBAR_MC_after * N_Jpsi_DDBAR_MC_after + \
            N_Jpsi_SSBAR_MC_after_relative_err * N_Jpsi_SSBAR_MC_after_relative_err * N_Jpsi_SSBAR_MC_after * N_Jpsi_SSBAR_MC_after + \
            N_Jpsi_CHARM_MC_after_relative_err * N_Jpsi_CHARM_MC_after_relative_err * N_Jpsi_CHARM_MC_after * N_Jpsi_CHARM_MC_after)) / (N_Jpsi_MC_after)
    else:
        N_Jpsi_MC_after_relative_err = 0
    Jpsi_MC_eff.append(efficiency( \
        N_Jpsi_SIGNAL_MC_after + \
        N_Jpsi_CHG_MC_after + \
        N_Jpsi_MIX_MC_after + \
        N_Jpsi_UUBAR_MC_after + \
        N_Jpsi_DDBAR_MC_after + \
        N_Jpsi_SSBAR_MC_after + \
        N_Jpsi_CHARM_MC_after, \
        N_Jpsi_SIGNAL_MC_before + \
        N_Jpsi_CHG_MC_before + \
        N_Jpsi_MIX_MC_before + \
        N_Jpsi_UUBAR_MC_before + \
        N_Jpsi_DDBAR_MC_before + \
        N_Jpsi_SSBAR_MC_before + \
        N_Jpsi_CHARM_MC_before))
    Jpsi_MC_eff_err.append((math.sqrt(
        N_Jpsi_MC_before_relative_err * N_Jpsi_MC_before_relative_err + N_Jpsi_MC_after_relative_err * N_Jpsi_MC_after_relative_err
        )) * \
        efficiency( \
        N_Jpsi_SIGNAL_MC_after + \
        N_Jpsi_CHG_MC_after + \
        N_Jpsi_MIX_MC_after + \
        N_Jpsi_UUBAR_MC_after + \
        N_Jpsi_DDBAR_MC_after + \
        N_Jpsi_SSBAR_MC_after + \
        N_Jpsi_CHARM_MC_after, \
        N_Jpsi_SIGNAL_MC_before + \
        N_Jpsi_CHG_MC_before + \
        N_Jpsi_MIX_MC_before + \
        N_Jpsi_UUBAR_MC_before + \
        N_Jpsi_DDBAR_MC_before + \
        N_Jpsi_SSBAR_MC_before + \
        N_Jpsi_CHARM_MC_before))

for N_embedded_data_before, N_embedded_data_after in zip(embedded_data_before, embedded_data_after):
    N_embedded_data_before_relative_err = relative_err_if_data(N_embedded_data_before)
    N_embedded_data_after_relative_err = relative_err_if_data(N_embedded_data_after)
    embedded_data_eff.append(efficiency(N_embedded_data_after, N_embedded_data_before))
    embedded_data_eff_err.append( \
        math.sqrt(N_embedded_data_before_relative_err*N_embedded_data_before_relative_err+N_embedded_data_after_relative_err*N_embedded_data_after_relative_err) * \
        efficiency(N_embedded_data_after, N_embedded_data_before)
        )

for N_embedded_CHG_MC_before, \
    N_embedded_MIX_MC_before, \
    N_embedded_UUBAR_MC_before, \
    N_embedded_DDBAR_MC_before, \
    N_embedded_SSBAR_MC_before, \
    N_embedded_CHARM_MC_before, \
    N_embedded_CHG_MC_after, \
    N_embedded_MIX_MC_after, \
    N_embedded_UUBAR_MC_after, \
    N_embedded_DDBAR_MC_after, \
    N_embedded_SSBAR_MC_after, \
    N_embedded_CHARM_MC_after \
    in zip(embedded_CHG_MC_before, \
           embedded_MIX_MC_before, \
           embedded_UUBAR_MC_before, \
           embedded_DDBAR_MC_before, \
           embedded_SSBAR_MC_before, \
           embedded_CHARM_MC_before, \
           embedded_CHG_MC_after, \
           embedded_MIX_MC_after, \
           embedded_UUBAR_MC_after, \
           embedded_DDBAR_MC_after, \
           embedded_SSBAR_MC_after, \
           embedded_CHARM_MC_after):
    N_embedded_CHG_MC_before_relative_err = relative_err_if_data(N_embedded_CHG_MC_before) * Scale_CHG_validation_MC15rd
    N_embedded_MIX_MC_before_relative_err = relative_err_if_data(N_embedded_MIX_MC_before) * Scale_MIX_validation_MC15rd
    N_embedded_UUBAR_MC_before_relative_err = relative_err_if_data(N_embedded_UUBAR_MC_before) * Scale_UUBAR_validation_MC15rd
    N_embedded_DDBAR_MC_before_relative_err = relative_err_if_data(N_embedded_DDBAR_MC_before) * Scale_DDBAR_validation_MC15rd
    N_embedded_SSBAR_MC_before_relative_err = relative_err_if_data(N_embedded_SSBAR_MC_before) * Scale_SSBAR_validation_MC15rd
    N_embedded_CHARM_MC_before_relative_err = relative_err_if_data(N_embedded_CHARM_MC_before) * Scale_CHARM_validation_MC15rd
    N_embedded_MC_before = \
                     N_embedded_CHG_MC_before + \
                     N_embedded_MIX_MC_before + \
                     N_embedded_UUBAR_MC_before + \
                     N_embedded_DDBAR_MC_before + \
                     N_embedded_SSBAR_MC_before + \
                     N_embedded_CHARM_MC_before
    if(N_embedded_MC_before != 0):
        N_embedded_MC_before_relative_err = (math.sqrt( \
            N_embedded_CHG_MC_before_relative_err * N_embedded_CHG_MC_before_relative_err * N_embedded_CHG_MC_before * N_embedded_CHG_MC_before + \
            N_embedded_MIX_MC_before_relative_err * N_embedded_MIX_MC_before_relative_err * N_embedded_MIX_MC_before * N_embedded_MIX_MC_before + \
            N_embedded_UUBAR_MC_before_relative_err * N_embedded_UUBAR_MC_before_relative_err * N_embedded_UUBAR_MC_before * N_embedded_UUBAR_MC_before + \
            N_embedded_DDBAR_MC_before_relative_err * N_embedded_DDBAR_MC_before_relative_err * N_embedded_DDBAR_MC_before * N_embedded_DDBAR_MC_before + \
            N_embedded_SSBAR_MC_before_relative_err * N_embedded_SSBAR_MC_before_relative_err * N_embedded_SSBAR_MC_before * N_embedded_SSBAR_MC_before + \
            N_embedded_CHARM_MC_before_relative_err * N_embedded_CHARM_MC_before_relative_err * N_embedded_CHARM_MC_before * N_embedded_CHARM_MC_before \
            )) / (N_embedded_MC_before)
    else:
         N_embedded_MC_before_relative_err = 0
    N_embedded_CHG_MC_after_relative_err = relative_err_if_data(N_embedded_CHG_MC_after) * Scale_CHG_validation_MC15rd
    N_embedded_MIX_MC_after_relative_err = relative_err_if_data(N_embedded_MIX_MC_after) * Scale_MIX_validation_MC15rd
    N_embedded_UUBAR_MC_after_relative_err = relative_err_if_data(N_embedded_UUBAR_MC_after) * Scale_UUBAR_validation_MC15rd
    N_embedded_DDBAR_MC_after_relative_err = relative_err_if_data(N_embedded_DDBAR_MC_after) * Scale_DDBAR_validation_MC15rd
    N_embedded_SSBAR_MC_after_relative_err = relative_err_if_data(N_embedded_SSBAR_MC_after) * Scale_SSBAR_validation_MC15rd
    N_embedded_CHARM_MC_after_relative_err = relative_err_if_data(N_embedded_CHARM_MC_after) * Scale_CHARM_validation_MC15rd
    N_embedded_MC_after = \
        N_embedded_CHG_MC_after + \
        N_embedded_MIX_MC_after + \
        N_embedded_UUBAR_MC_after + \
        N_embedded_DDBAR_MC_after + \
        N_embedded_SSBAR_MC_after + \
        N_embedded_CHARM_MC_after
    if(N_embedded_MC_after != 0):
        N_embedded_MC_after_relative_err = (math.sqrt( \
            N_embedded_CHG_MC_after_relative_err * N_embedded_CHG_MC_after_relative_err * N_embedded_CHG_MC_after * N_embedded_CHG_MC_after + \
            N_embedded_MIX_MC_after_relative_err * N_embedded_MIX_MC_after_relative_err * N_embedded_MIX_MC_after * N_embedded_MIX_MC_after + \
            N_embedded_UUBAR_MC_after_relative_err * N_embedded_UUBAR_MC_after_relative_err * N_embedded_UUBAR_MC_after * N_embedded_UUBAR_MC_after + \
            N_embedded_DDBAR_MC_after_relative_err * N_embedded_DDBAR_MC_after_relative_err * N_embedded_DDBAR_MC_after * N_embedded_DDBAR_MC_after + \
            N_embedded_SSBAR_MC_after_relative_err * N_embedded_SSBAR_MC_after_relative_err * N_embedded_SSBAR_MC_after * N_embedded_SSBAR_MC_after + \
            N_embedded_CHARM_MC_after_relative_err * N_embedded_CHARM_MC_after_relative_err * N_embedded_CHARM_MC_after * N_embedded_CHARM_MC_after)) / (N_embedded_MC_after)
    else:
        N_embedded_MC_after_relative_err = 0
    embedded_MC_eff.append(efficiency( \
        N_embedded_CHG_MC_after + \
        N_embedded_MIX_MC_after + \
        N_embedded_UUBAR_MC_after + \
        N_embedded_DDBAR_MC_after + \
        N_embedded_SSBAR_MC_after + \
        N_embedded_CHARM_MC_after, \
        N_embedded_CHG_MC_before + \
        N_embedded_MIX_MC_before + \
        N_embedded_UUBAR_MC_before + \
        N_embedded_DDBAR_MC_before + \
        N_embedded_SSBAR_MC_before + \
        N_embedded_CHARM_MC_before))
    embedded_MC_eff_err.append((math.sqrt(
        N_embedded_MC_before_relative_err * N_embedded_MC_before_relative_err + N_embedded_MC_after_relative_err * N_embedded_MC_after_relative_err
        )) * \
        efficiency( \
        N_embedded_CHG_MC_after + \
        N_embedded_MIX_MC_after + \
        N_embedded_UUBAR_MC_after + \
        N_embedded_DDBAR_MC_after + \
        N_embedded_SSBAR_MC_after + \
        N_embedded_CHARM_MC_after, \
        N_embedded_CHG_MC_before + \
        N_embedded_MIX_MC_before + \
        N_embedded_UUBAR_MC_before + \
        N_embedded_DDBAR_MC_before + \
        N_embedded_SSBAR_MC_before + \
        N_embedded_CHARM_MC_before))

plot_and_save_data(x, SIGNAL_nominal_MC_eff, Jpsi_data_eff, Jpsi_data_eff_err, Jpsi_MC_eff, Jpsi_MC_eff_err, embedded_data_eff, embedded_data_eff_err, embedded_MC_eff, embedded_MC_eff_err, "FBDT_eff_compare.png", "")


#plot_and_save_data(x, y1_1, yerr1_1, "FBDT_eff.png", r'$B \rightarrow X_{s} J/\psi$', r'')
