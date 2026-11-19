"""Plot combined hadronic FEI calibration factors for P_FEI > 0.001.

Values and total uncertainties are transcribed at the published precision from
Table I (page 8) of BELLE2-NOTE-PH-2024-036 (Version 3.0, July 26, 2024).
Both charged and neutral tags use the 0.001 threshold. No combined "All"
factor is provided in that table, so no inclusive row is drawn.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 15,
    "axes.labelsize": 15,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "pdf.fonttype": 42,
})

CHARGED_LABELS = [
    r"$\bar{D}^0\pi^+$", r"$\bar{D}^0\pi^+\pi^0$",
    r"$\bar{D}^0\pi^+\pi^+\pi^-$", r"$\bar{D}^0\pi^+\pi^+\pi^-\pi^0$",
    r"$\bar{D}^{*0}\pi^+$", r"$\bar{D}^{*0}\pi^+\pi^0$",
    r"$\bar{D}^{*0}\pi^+\pi^+\pi^-$", r"$\bar{D}^{*0}\pi^+\pi^+\pi^-\pi^0$",
    r"$D^-\pi^+\pi^+$", r"$D^-\pi^+\pi^+\pi^0$",
    r"$\Lambda_c^- p\pi^+\pi^-\pi^+$", "rest",
]
CHARGED_VALUES = [1.06, 0.82, 0.78, 0.67, 1.00, 1.01, 0.84, 0.68, 0.88, 0.64, 0.38, 0.77]
CHARGED_ERRORS = [0.03, 0.06, 0.07, 0.16, 0.03, 0.03, 0.03, 0.08, 0.20, 0.09, 0.12, 0.04]

NEUTRAL_LABELS = [
    r"$D^-\pi^+$", r"$D^-\pi^+\pi^0$",
    r"$D^-\pi^+\pi^+\pi^-$", r"$D^-\pi^+\pi^+\pi^-\pi^0$",
    r"$\bar{D}^0\pi^+\pi^-$", r"$D^{*-}\pi^+$",
    r"$D^{*-}\pi^+\pi^0$", r"$D^{*-}\pi^+\pi^+\pi^-$",
    r"$D^{*-}\pi^+\pi^+\pi^-\pi^0$",
    r"$\Lambda_c^- p\pi^+\pi^-$", "rest",
]
NEUTRAL_VALUES = [1.12, 1.01, 0.88, 0.89, 1.06, 1.02, 0.68, 0.71, 0.72, 0.60, 0.91]
NEUTRAL_ERRORS = [0.04, 0.03, 0.08, 0.03, 0.03, 0.03, 0.02, 0.02, 0.09, 0.02, 0.03]


def make_plot(labels, values, errors, tag):
    """Keep the original horizontal-errorbar style, with one combined series."""
    if not len(labels) == len(values) == len(errors):
        raise ValueError("Labels, factors, and uncertainties must have equal lengths.")

    fig, ax = plt.subplots(figsize=(8, 6))
    y_positions = np.arange(len(labels))
    ax.errorbar(
        values[::-1], y_positions, xerr=errors[::-1],
        fmt="o", color="black",
    )
    ax.set_yticks(y_positions, labels=labels[::-1])
    ax.set_xlabel("Calibration Factor")
    ax.set_xlim(0, 1.4)
    ax.set_title(r"$\mathcal{P}_{\mathrm{FEI}} > 0.001, " + tag + r"$")
    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--no-show", action="store_true", help="Save plots without opening windows.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plots = [
        ("Bp", CHARGED_LABELS, CHARGED_VALUES, CHARGED_ERRORS, "B^+"),
        ("B0", NEUTRAL_LABELS, NEUTRAL_VALUES, NEUTRAL_ERRORS, "B^0"),
    ]
    for suffix, labels, values, errors, tag in plots:
        fig = make_plot(labels, values, errors, tag)
        for extension in ("pdf", "png"):
            output = args.output_dir / f"FEI_cal_combined_{suffix}.{extension}"
            fig.savefig(output, dpi=300, facecolor="white")
            print(f"Saved: {output.resolve()}")
        if args.no_show:
            plt.close(fig)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
