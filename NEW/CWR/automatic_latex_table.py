import os
import re
import math

def extract_minos_for_param(filepath, param_name):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if f"fit result for {param_name}" in line:
                    if i + 1 < len(lines) and "MINOS error:" in lines[i + 1]:
                        match = re.search(r"MINOS error:\s*([+-]?\d*\.\d+)\s+([+-]?\d*\.\d+)", lines[i + 1])
                        if match:
                            return float(match.group(1)), float(match.group(2))
    except FileNotFoundError:
        pass
    return None, None

def compute_uncertainty_contribution(fix_none, fix_x):
    if fix_none[0] is None or fix_x[0] is None:
        return None
    uncer = math.sqrt(abs(((abs(fix_none[0]) + abs(fix_none[1]))/2.0)**2 - ((abs(fix_x[0]) + abs(fix_x[1]))/2.0)**2))
    return uncer

# Mapping of LaTeX Belle note label, paper label, type, size, and directory name
systematics = [
    ("Background normalization", "Background normalization", "\\texttt{OverallSys}", "20\\%", "fix_BKGNorm"),
    ("MC statistics", "MC statistics", "\\texttt{StatError}", "$O(1\\%)$", "fix_MCstat"),
    ("Photon multiplicity correction", "Photon multiplicity correction", "\\texttt{HistoSys}", "$O(1\\%)$", "fix_photon_multiplicity"),
    ("BR of main B meson decays", "Branching ratio of major $B$ meson decay", "\\texttt{HistoSys}", "$O(1 \\text{-} 100\\%)$", "fix_BBBR"),
    ("$q\\bar{q}$ shape", "${\\mathcal{O}}$ distribution shape for $q\\bar{q}$ background", "\\texttt{HistoSys}", "100\% of weight", "fix_BDTc"),
    ("$q\\bar{q}$ efficiency", "$q\\bar{q}$ background efficiency", "\\texttt{OverallSys}", "$O(10\\%)$", "fix_qqbar_CAL"),
    ("FastBDT efficiency", "$\\mathcal{O}$ selection efficiency", "\\texttt{OverallSys}", "$O(1 \\text{-} 10\\%)$", "fix_FBDT_efficiency"),
    ("$K^{*}-X_{s}$ transition", "Non-resonant $X_{s}\\nu\\bar{\\nu}$ generation point", "\\tiny \\texttt{HistoSys} and \\texttt{OverallSys}", "9.1\\%", "fix_Transition"),
    ("Fragmentation", "Fragmentation", "\\texttt{HistoSys}", "$O(1 \\text{-} 100\\%)$", "fix_Fragmentation"),
    ("Fraction of decay modes", "Fraction of $K\\nu\\bar{\\nu}$, $K^{*}\\nu\\bar{\\nu}$ decays", "\\texttt{HistoSys}", "$O(5\\%)$", "fix_fraction"),
    ("Another FEI calibration factor", "Tagging efficiency from $B \\to K^{(*)}J/\\psi$ decay", "\\texttt{HistoSys}", "$O(20 \\text{-} 100\\%)$", "fix_NEWFEICAL"),
    ("FEI calibration for $B\\bar{B}$", "Tagging efficiency", "\\texttt{HistoSys}", "$O(1 \\text{-} 30\\%)$", "fix_FEI"),
    ("$\\pi^{0}$ reconstruction", "$\\pi^{0}$ reconstruction efficiency", "\\texttt{HistoSys}", "$O(5\\%)$", "fix_pi0"),
    ("$B\\rightarrow K$ form factor", "$B\\to K$ form factor", "\\texttt{HistoSys}", "$O(1\\%)$", "fix_Kff"),
    ("$K_{S}^{0}$ reconstruction", "$K_{S}^{0}$ reconstruction efficiency", "\\texttt{HistoSys}", "$O(5\\%)$", "fix_KS0"),
    ("$B\\rightarrow K^{*}$ form factor", "$B\\to K^{*}$ form factor", "\\texttt{HistoSys}", "$O(10\\%)$", "fix_Kstarff"),
    ("BR of $B \\rightarrow X_{s} K_{L}^{0} K_{L}^{0}$", "Branching ratio for $B \\to X_{s} K_{L}^{0} K_{L}^{0}$", "\\texttt{HistoSys}", "100\\%", "fix_BRBtoXKLKL"),
    ("BB counting", "The number of $B\\bar{B}$ pair", "\\texttt{OverallSys}", "1.45\\%", "fix_BBcounting"),
    ("Pion ID", "Charged pion identification", "\\texttt{HistoSys}", "$O(1\\%)$", "fix_PID"),
    ("Kaon ID", "Charged kaon identification", "\\texttt{HistoSys}", "$O(1\\%)$", "fix_KID"),
    ("BR of $D \\rightarrow K_{L}^{0} X$ from B meson", "Branching ratio for $B \\to X_{s} K_{L}^{0} K_{L}^{0}$", "\\texttt{HistoSys}", "7.7\\%", "fix_BtoDtoXKL"),
    ("BR of $B\\rightarrow X_{s} n\\bar{n}$ decay", "Branching ratio for $B \\to X_{s} n\\bar{n}$", "\\texttt{HistoSys}", "100\\%", "fix_BRXnn"),
    ("Fermi motion momentum", "Fermi motion momentum", "\\texttt{HistoSys}", "$O(8\\%)$", "fix_pf"),
    ("b-quark mass", "b-quark mass parameter", "\\texttt{HistoSys}", "3.2\\%", "fix_mb"),
    ("Tracking efficiency", "Tracking efficiency", "\\texttt{HistoSys}", "0.24\\%", "fix_track"),
    ("$K_{L}^{0}$ ECL efficiency", "$K_{L}^{0}$ efficiency in the ECL", "\\texttt{HistoSys}", "17\\%", "fix_EffECLKL"),
    ("fitter bias", "Fitting bias", "$\\texttt{OverallSys}$", "$O(5\\%)$", "fix_fitter"),
]

IsItFullUnblindedDir = True

if (IsItFullUnblindedDir == False):
    # ========================================= Belle note one mu ========================================= #
    print("Belle note one mu:")

    # Extract base errors
    param = "mu"
    fix_none = extract_minos_for_param("./one_mu/fix_none/log/FIX_1.log", param)
    fix_all = extract_minos_for_param("./one_mu/fix_all/log/FIX_1.log", param)

    # Generate LaTeX table rows
    unsorted_rows = []
    sorted_rows = []
    latex_rows = []

    for label_note, label_paper, systype, syssize, dir in systematics:
        fix_dir = extract_minos_for_param(f"./one_mu/{dir}/log/FIX_1.log", param)
        contrib = compute_uncertainty_contribution(fix_none, fix_dir)
        if contrib is not None:
            unsorted_rows.append((label_note, systype, contrib))

    # lets sort
    sorted_rows = sorted(unsorted_rows, key=lambda x: x[2], reverse=True)

    # convert into latex rows
    for label_note, systype, contrib in sorted_rows:
        row = f"{label_note} & ${systype}$ & {contrib:.2f} \\\\"
        latex_rows.append(row)

    # Append statistical uncertainty from fix_all
    if fix_all[0] is not None and fix_all[1] is not None:
        value = (abs(fix_all[0]) + abs(fix_all[1]))/2.0
        stat_row = f"\\hline\nstatistical uncertainty & & {value:.2f}\\\\"
        latex_rows.append(stat_row)

    # Print LaTeX table
    print("\n".join(latex_rows))
    print("\n")
    # ===================================================================================================== #

    # ========================================= Belle note mu ========================================= #
    print("Belle note free_all_mu:")

    # Extract base errors
    param_mu = "mu"
    param_MXs1 = "mu_MXs1"
    param_MXs2 = "mu_MXs2"
    param_MXs3 = "mu_MXs3"

    fix_none_mu = extract_minos_for_param("./free_all_mu/fix_none/log/FIX_1.log", param_mu)
    fix_all_mu = extract_minos_for_param("./free_all_mu/fix_all/log/FIX_1.log", param_mu)
    fix_none_MXs1 = extract_minos_for_param("./free_all/fix_none/log/FIX_1.log", param_MXs1)
    fix_all_MXs1 = extract_minos_for_param("./free_all/fix_all/log/FIX_1.log", param_MXs1)
    fix_none_MXs2 = extract_minos_for_param("./free_all/fix_none/log/FIX_1.log", param_MXs2)
    fix_all_MXs2 = extract_minos_for_param("./free_all/fix_all/log/FIX_1.log", param_MXs2)
    fix_none_MXs3 = extract_minos_for_param("./free_all/fix_none/log/FIX_1.log", param_MXs3)
    fix_all_MXs3 = extract_minos_for_param("./free_all/fix_all/log/FIX_1.log", param_MXs3)

    # Generate LaTeX table rows
    unsorted_rows = []
    sorted_rows = []
    latex_rows = []

    for label_note, label_paper, systype, syssize, dir in systematics:
        fix_dir_mu = extract_minos_for_param(f"./free_all_mu/{dir}/log/FIX_1.log", param_mu)
        contrib_mu = compute_uncertainty_contribution(fix_none_mu, fix_dir_mu)
        fix_dir_MXs1 = extract_minos_for_param(f"./free_all/{dir}/log/FIX_1.log", param_MXs1)
        contrib_MXs1 = compute_uncertainty_contribution(fix_none_MXs1, fix_dir_MXs1)
        fix_dir_MXs2 = extract_minos_for_param(f"./free_all/{dir}/log/FIX_1.log", param_MXs2)
        contrib_MXs2 = compute_uncertainty_contribution(fix_none_MXs2, fix_dir_MXs2)
        fix_dir_MXs3 = extract_minos_for_param(f"./free_all/{dir}/log/FIX_1.log", param_MXs3)
        contrib_MXs3 = compute_uncertainty_contribution(fix_none_MXs3, fix_dir_MXs3)
        
        if contrib_MXs1 is not None:
            unsorted_rows.append((label_note, systype, contrib_MXs1, contrib_MXs2, contrib_MXs3))

    # lets sort
    sorted_rows = sorted(unsorted_rows, key=lambda x: x[4], reverse=True)

    # convert into latex rows
    for label_note, systype, contrib_MXs1, contrib_MXs2, contrib_MXs3 in sorted_rows:
        row = f"{label_note} & {systype} & {contrib_MXs1:.2f} & {contrib_MXs2:.2f} & {contrib_MXs3:.2f} \\\\"
        latex_rows.append(row)

    # Append statistical uncertainty from fix_all
    if fix_all_MXs1[0] is not None and fix_all_MXs1[1] is not None and fix_all_MXs2[0] is not None and fix_all_MXs2[1] is not None and fix_all_MXs3[0] is not None and fix_all_MXs3[1] is not None:
        value_MXs1 = (abs(fix_all_MXs1[0]) + abs(fix_all_MXs1[1]))/2.0
        value_MXs2 = (abs(fix_all_MXs2[0]) + abs(fix_all_MXs2[1]))/2.0
        value_MXs3 = (abs(fix_all_MXs3[0]) + abs(fix_all_MXs3[1]))/2.0
        stat_row = f"\\hline\nStatistical uncertainty & & {value_MXs1:.2f} & {value_MXs2:.2f} & {value_MXs3:.2f} \\\\"
        latex_rows.append(stat_row)

    # Print LaTeX table
    print("\n".join(latex_rows))
    print("\n")
    # ===================================================================================================== #
    
else:# in full_unblind directory
    # ========================================= paper mu in main text ========================================= #
    print("paper free_all_mu in main text:")

    # Extract base errors
    param = "mu"
    BR_mu = 0.000029
    Scale = 100000.0

    fix_none = extract_minos_for_param("./free_all_fitter_bias_mu/fix_none/log/FIX_1.log", param)
    fix_all = extract_minos_for_param("./free_all_fitter_bias_mu/fix_all/log/FIX_1.log", param)
    fix_minor = extract_minos_for_param("./free_all_fitter_bias_mu/fix_minor/log/FIX_1.log", param)

    # Generate LaTeX table rows
    unsorted_rows = []
    sorted_rows = []
    latex_rows = []

    for label_note, label_paper, systype, syssize, dir in systematics:
        fix_dir = extract_minos_for_param(f"./free_all_fitter_bias_mu/{dir}/log/FIX_1.log", param)
        contrib = compute_uncertainty_contribution(fix_none, fix_dir)
        if contrib is not None:
            value = Scale * BR_mu * contrib
            unsorted_rows.append((label_paper, value))

    # lets sort
    sorted_rows = sorted(unsorted_rows, key=lambda x: x[1], reverse=True)

    # convert into latex rows (only 7 element)
    for label_paper, impact in sorted_rows[:7]:
        row = f"{label_paper} & {impact:.1f}\\\\"
        latex_rows.append(row)

    # minor uncertainty
    if fix_minor[0] is not None and fix_minor[1] is not None:
        contrib = compute_uncertainty_contribution(fix_none, fix_minor)
        if contrib is not None:
            value = Scale * BR_mu * contrib
            row = f"Other subdominant contributions & {value:.1f}\\\\"
            latex_rows.append(row)

    # Append total systematic uncertainty from fix_all
    if fix_all[0] is not None and fix_all[1] is not None:
        contrib = compute_uncertainty_contribution(fix_none, fix_all)
        if contrib is not None:
            value = Scale * BR_mu * contrib
            row = f"\\hline\nTotal systematic sources & {value:.1f}\\\\"
            latex_rows.append(row)

    # Print LaTeX table
    print("\n".join(latex_rows))
    print("\n")
    # ===================================================================================================== #

    # ========================================= Belle note one mu ========================================= #
    print("Belle note one mu:")

    # Extract base errors
    param = "mu"
    fix_none = extract_minos_for_param("./one_mu_fitter_bias/fix_none/log/FIX_1.log", param)
    fix_all = extract_minos_for_param("./one_mu_fitter_bias/fix_all/log/FIX_1.log", param)

    # Generate LaTeX table rows
    unsorted_rows = []
    sorted_rows = []
    latex_rows = []

    for label_note, label_paper, systype, syssize, dir in systematics:
        fix_dir = extract_minos_for_param(f"./one_mu_fitter_bias/{dir}/log/FIX_1.log", param)
        contrib = compute_uncertainty_contribution(fix_none, fix_dir)
        if contrib is not None:
            unsorted_rows.append((label_note, systype, contrib))

    # lets sort
    sorted_rows = sorted(unsorted_rows, key=lambda x: x[2], reverse=True)

    # convert into latex rows
    for label_note, systype, contrib in sorted_rows:
        row = f"{label_note} & ${systype}$ & {contrib:.2f} \\\\"
        latex_rows.append(row)

    # Append statistical uncertainty from fix_all
    if fix_all[0] is not None and fix_all[1] is not None:
        value = (abs(fix_all[0]) + abs(fix_all[1]))/2.0
        stat_row = f"\\hline\nstatistical uncertainty & & {value:.2f}\\\\"
        latex_rows.append(stat_row)

    # Print LaTeX table
    print("\n".join(latex_rows))
    print("\n")
    # ===================================================================================================== #

    # ========================================= paper mu in supplemental material ========================================= #
    print("paper free_all_mu in supplemental material:")

    # Extract base errors
    param_mu = "mu"
    param_MXs1 = "mu_MXs1"
    param_MXs2 = "mu_MXs2"
    param_MXs3 = "mu_MXs3"

    BR_mu = 0.000029
    BR_MXs1 = 0.0000048502
    BR_MXs2 = 0.0000085007
    BR_MXs3 = 0.0000156683
    Scale = 100000.0

    fix_none_mu = extract_minos_for_param("./free_all_fitter_bias_mu/fix_none/log/FIX_1.log", param_mu)
    fix_all_mu = extract_minos_for_param("./free_all_fitter_bias_mu/fix_all/log/FIX_1.log", param_mu)
    fix_none_MXs1 = extract_minos_for_param("./free_all_fitter_bias/fix_none/log/FIX_1.log", param_MXs1)
    fix_all_MXs1 = extract_minos_for_param("./free_all_fitter_bias/fix_all/log/FIX_1.log", param_MXs1)
    fix_none_MXs2 = extract_minos_for_param("./free_all_fitter_bias/fix_none/log/FIX_1.log", param_MXs2)
    fix_all_MXs2 = extract_minos_for_param("./free_all_fitter_bias/fix_all/log/FIX_1.log", param_MXs2)
    fix_none_MXs3 = extract_minos_for_param("./free_all_fitter_bias/fix_none/log/FIX_1.log", param_MXs3)
    fix_all_MXs3 = extract_minos_for_param("./free_all_fitter_bias/fix_all/log/FIX_1.log", param_MXs3)

    # Generate LaTeX table rows
    unsorted_rows = []
    sorted_rows = []
    latex_rows = []

    for label_note, label_paper, systype, syssize, dir in systematics:
        fix_dir_mu = extract_minos_for_param(f"./free_all_fitter_bias_mu/{dir}/log/FIX_1.log", param_mu)
        contrib_mu = compute_uncertainty_contribution(fix_none_mu, fix_dir_mu)
        fix_dir_MXs1 = extract_minos_for_param(f"./free_all_fitter_bias/{dir}/log/FIX_1.log", param_MXs1)
        contrib_MXs1 = compute_uncertainty_contribution(fix_none_MXs1, fix_dir_MXs1)
        fix_dir_MXs2 = extract_minos_for_param(f"./free_all_fitter_bias/{dir}/log/FIX_1.log", param_MXs2)
        contrib_MXs2 = compute_uncertainty_contribution(fix_none_MXs2, fix_dir_MXs2)
        fix_dir_MXs3 = extract_minos_for_param(f"./free_all_fitter_bias/{dir}/log/FIX_1.log", param_MXs3)
        contrib_MXs3 = compute_uncertainty_contribution(fix_none_MXs3, fix_dir_MXs3)
        
        if contrib is not None:
            value_mu = Scale * BR_mu * contrib_mu
            value_MXs1 = Scale * BR_MXs1 * contrib_MXs1
            value_MXs2 = Scale * BR_MXs2 * contrib_MXs2
            value_MXs3 = Scale * BR_MXs3 * contrib_MXs3
        
            unsorted_rows.append((label_paper, syssize, value_mu, value_MXs1, value_MXs2, value_MXs3))

    # lets sort
    sorted_rows = sorted(unsorted_rows, key=lambda x: x[2], reverse=True)

    # convert into latex rows
    for label_paper, syssize, impact_mu, impact_MXs1, impact_MXs2, impact_MXs3 in sorted_rows:
        row = f"{label_paper} & {syssize} & {impact_mu:.1f} & {impact_MXs1:.1f} & {impact_MXs2:.1f} & {impact_MXs3:.1f} \\\\"
        latex_rows.append(row)

    # Append all systematic uncertainty from fix_all
    if fix_all_mu[0] is not None and fix_all_mu[1] is not None and fix_all_MXs1[0] is not None and fix_all_MXs1[1] is not None and fix_all_MXs2[0] is not None and fix_all_MXs2[1] is not None and fix_all_MXs3[0] is not None and fix_all_MXs3[1] is not None:
        value_mu = Scale * BR_mu * compute_uncertainty_contribution(fix_none_mu, fix_all_mu)
        value_MXs1 = Scale * BR_MXs1 * compute_uncertainty_contribution(fix_none_MXs1, fix_all_MXs1)
        value_MXs2 = Scale * BR_MXs2 * compute_uncertainty_contribution(fix_none_MXs2, fix_all_MXs2)
        value_MXs3 = Scale * BR_MXs3 * compute_uncertainty_contribution(fix_none_MXs3, fix_all_MXs3)
        syst_row = f"\\hline\nTotal systematic sources & & {value_mu:.1f} & {value_MXs1:.1f} & {value_MXs2:.1f} & {value_MXs3:.1f} \\\\"
        latex_rows.append(syst_row)

    # Append statistical uncertainty from fix_all
    if fix_all_mu[0] is not None and fix_all_mu[1] is not None and fix_all_MXs1[0] is not None and fix_all_MXs1[1] is not None and fix_all_MXs2[0] is not None and fix_all_MXs2[1] is not None and fix_all_MXs3[0] is not None and fix_all_MXs3[1] is not None:
        value_mu = Scale * BR_mu * (abs(fix_all_mu[0]) + abs(fix_all_mu[1]))/2.0
        value_MXs1 = Scale * BR_MXs1 * (abs(fix_all_MXs1[0]) + abs(fix_all_MXs1[1]))/2.0
        value_MXs2 = Scale * BR_MXs2 * (abs(fix_all_MXs2[0]) + abs(fix_all_MXs2[1]))/2.0
        value_MXs3 = Scale * BR_MXs3 * (abs(fix_all_MXs3[0]) + abs(fix_all_MXs3[1]))/2.0
        stat_row = f"\\hline\nStatistical uncertainty & & {value_mu:.1f} & {value_MXs1:.1f} & {value_MXs2:.1f} & {value_MXs3:.1f} \\\\"
        latex_rows.append(stat_row)

    # Print LaTeX table
    print("\n".join(latex_rows))
    print("\n")
    # ===================================================================================================== #

    # ========================================= Belle note mu ========================================= #
    print("Belle note free_all_mu:")

    # Extract base errors
    param_mu = "mu"
    param_MXs1 = "mu_MXs1"
    param_MXs2 = "mu_MXs2"
    param_MXs3 = "mu_MXs3"

    fix_none_mu = extract_minos_for_param("./free_all_fitter_bias_mu/fix_none/log/FIX_1.log", param_mu)
    fix_all_mu = extract_minos_for_param("./free_all_fitter_bias_mu/fix_all/log/FIX_1.log", param_mu)
    fix_none_MXs1 = extract_minos_for_param("./free_all_fitter_bias/fix_none/log/FIX_1.log", param_MXs1)
    fix_all_MXs1 = extract_minos_for_param("./free_all_fitter_bias/fix_all/log/FIX_1.log", param_MXs1)
    fix_none_MXs2 = extract_minos_for_param("./free_all_fitter_bias/fix_none/log/FIX_1.log", param_MXs2)
    fix_all_MXs2 = extract_minos_for_param("./free_all_fitter_bias/fix_all/log/FIX_1.log", param_MXs2)
    fix_none_MXs3 = extract_minos_for_param("./free_all_fitter_bias/fix_none/log/FIX_1.log", param_MXs3)
    fix_all_MXs3 = extract_minos_for_param("./free_all_fitter_bias/fix_all/log/FIX_1.log", param_MXs3)

    # Generate LaTeX table rows
    unsorted_rows = []
    sorted_rows = []
    latex_rows = []

    for label_note, label_paper, systype, syssize, dir in systematics:
        fix_dir_mu = extract_minos_for_param(f"./free_all_fitter_bias_mu/{dir}/log/FIX_1.log", param_mu)
        contrib_mu = compute_uncertainty_contribution(fix_none_mu, fix_dir_mu)
        fix_dir_MXs1 = extract_minos_for_param(f"./free_all_fitter_bias/{dir}/log/FIX_1.log", param_MXs1)
        contrib_MXs1 = compute_uncertainty_contribution(fix_none_MXs1, fix_dir_MXs1)
        fix_dir_MXs2 = extract_minos_for_param(f"./free_all_fitter_bias/{dir}/log/FIX_1.log", param_MXs2)
        contrib_MXs2 = compute_uncertainty_contribution(fix_none_MXs2, fix_dir_MXs2)
        fix_dir_MXs3 = extract_minos_for_param(f"./free_all_fitter_bias/{dir}/log/FIX_1.log", param_MXs3)
        contrib_MXs3 = compute_uncertainty_contribution(fix_none_MXs3, fix_dir_MXs3)
        
        if contrib is not None:
            unsorted_rows.append((label_note, systype, contrib_MXs1, contrib_MXs2, contrib_MXs3))

    # lets sort
    sorted_rows = sorted(unsorted_rows, key=lambda x: x[4], reverse=True)

    # convert into latex rows
    for label_note, systype, contrib_MXs1, contrib_MXs2, contrib_MXs3 in sorted_rows:
        row = f"{label_note} & {systype} & {contrib_MXs1:.2f} & {contrib_MXs2:.2f} & {contrib_MXs3:.2f} \\\\"
        latex_rows.append(row)

    # Append statistical uncertainty from fix_all
    if fix_all_MXs1[0] is not None and fix_all_MXs1[1] is not None and fix_all_MXs2[0] is not None and fix_all_MXs2[1] is not None and fix_all_MXs3[0] is not None and fix_all_MXs3[1] is not None:
        value_mu = (abs(fix_all_mu[0]) + abs(fix_all_mu[1]))/2.0
        value_MXs1 = (abs(fix_all_MXs1[0]) + abs(fix_all_MXs1[1]))/2.0
        value_MXs2 = (abs(fix_all_MXs2[0]) + abs(fix_all_MXs2[1]))/2.0
        value_MXs3 = (abs(fix_all_MXs3[0]) + abs(fix_all_MXs3[1]))/2.0
        stat_row = f"\\hline\nStatistical uncertainty & & {value_MXs1:.2f} & {value_MXs2:.2f} & {value_MXs3:.2f} \\\\"
        latex_rows.append(stat_row)

    # Print LaTeX table
    print("\n".join(latex_rows))
    print("\n")
    # ===================================================================================================== #

