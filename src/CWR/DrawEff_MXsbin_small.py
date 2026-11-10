"""Plot the four efficiency comparisons used in thesis Figs. M.1--M.4.

Run without arguments to save all four modes as PNG and vector PDF files
beside this script. Use --mode for one mode or --output-dir for another
existing directory. Input values and NaNs are retained from the original
DrawEff_MXsbin_small.py; a NaN omits the point rather than plotting zero.

With an 8-inch canvas inserted at 0.7\\textwidth in the thesis, 16-point
labels become approximately 8.8 points in the compiled document.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


X_BOUNDARIES = np.array([1.15, 1.5, 2.0, 2.4, 3.0])

DATASETS = {
    "Kcpicpic": {  # Fig. M.1
        "pythia": [0.001650, 0.001382, 0.000303, 0.000032],
        "pythia_error": [0.000042, 0.000044, 0.000044, 0.000032],
        "k1": [0.001714, 0.001413, 0.000154, np.nan],
        "k1_error": [0.000076, 0.000215, 0.000154, np.nan],
        "description": r"$B^{+}\rightarrow K^{+}\pi^{-}\pi^{+}\nu\bar{\nu}$",
    },
    "K0picpi0": {  # Fig. M.2
        "pythia": [0.000522, 0.000495, 0.000283, 0.000052],
        "pythia_error": [0.000022, 0.000023, 0.000038, 0.000037],
        "k1": [0.000618, 0.000484, np.nan, np.nan],
        "k1_error": [0.000040, 0.000099, np.nan, np.nan],
        "description": r"$B^{+}\rightarrow K^{0}\pi^{+}\pi^{0}\nu\bar{\nu}$",
    },
    "K0picpic": {  # Fig. M.3
        "pythia": [0.000225, 0.000203, 0.000057, np.nan],
        "pythia_error": [0.000016, 0.000017, 0.000019, np.nan],
        "k1": [0.000295, 0.000166, np.nan, np.nan],
        "k1_error": [0.000031, 0.000074, np.nan, np.nan],
        "description": r"$B^{0}\rightarrow K^{0}\pi^{+}\pi^{-}\nu\bar{\nu}$",
    },
    "Kcpicpi0": {  # Fig. M.4
        "pythia": [0.001389, 0.001219, 0.000493, 0.000513],
        "pythia_error": [0.000035, 0.000035, 0.000050, 0.000115],
        "k1": [0.001439, 0.001299, 0.000455, np.nan],
        "k1_error": [0.000060, 0.000162, 0.000204, np.nan],
        "description": r"$B^{0}\rightarrow K^{+}\pi^{-}\pi^{0}\nu\bar{\nu}$",
    },
}

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 16,
    "mathtext.fontset": "dejavusans",
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def make_figure(mode):
    """Return one figure without writing files."""
    data = DATASETS[mode]
    x_centers = (X_BOUNDARIES[:-1] + X_BOUNDARIES[1:]) / 2
    x_errors = np.diff(X_BOUNDARIES) / 2

    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(
            x_centers, data["pythia"],
            xerr=x_errors, yerr=data["pythia_error"],
            fmt="o", color="tab:blue", capsize=5,
            label=r"$B\rightarrow X_s\nu\bar{\nu}$ MC with PYTHIA",
        )
        ax.errorbar(
            x_centers, data["k1"],
            xerr=x_errors, yerr=data["k1_error"],
            fmt="s", color="tab:orange", capsize=5,
            label=r"$B\rightarrow K_1(1270)\nu\bar{\nu}$ MC",
        )

        for boundary in X_BOUNDARIES:
            ax.axvline(
                boundary, color="gray", linestyle="--", linewidth=0.5,
            )

        ax.set_xlabel(r"$M_{X_s}$ [GeV/$c^2$]")
        ax.set_ylabel("efficiency")
        ax.set_xticks(np.arange(1.25, 3.01, 0.25))
        ax.legend(
            loc="upper right", borderpad=0.3,
            handlelength=1.0, handletextpad=0.4, labelspacing=0.3,
        )
        ax.text(
            0.97, 0.74, data["description"],
            transform=ax.transAxes, fontsize=18,
            ha="right", va="top",
        )
        fig.tight_layout(pad=0.6)

    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode", choices=["all", *DATASETS], default="all",
        help="Mode to plot (default: all four modes).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent,
        help="Existing output directory (default: this script's directory).",
    )
    args = parser.parse_args()
    if not args.output_dir.is_dir():
        parser.error(f"Output directory does not exist: {args.output_dir}")

    modes = DATASETS if args.mode == "all" else [args.mode]
    for mode in modes:
        fig = make_figure(mode)
        try:
            with plt.rc_context(STYLE):
                for extension in ("png", "pdf"):
                    destination = args.output_dir / f"{mode}.{extension}"
                    # Keep the fixed canvas size so all panels have the same scale.
                    fig.savefig(destination, dpi=600, facecolor="white")
                    print(destination)
        finally:
            plt.close(fig)


if __name__ == "__main__":
    main()
