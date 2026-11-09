"""Replot the kaon PID correction table shown in thesis Fig. K.1.

The values and uncertainties are transcribed from the displayed labels in
phd_thesis/picture/systematic/KaonEff.png. These are rounded to two decimals;
they are NOT the original calibration data. In particular, a displayed
uncertainty of 0.00 must not be interpreted as an exact zero uncertainty.
Use these arrays only to reproduce the figure, not for PID corrections.

The original plasma colour scale (0 to 1), equal-size bin cells, downward
momentum axis and blue hatching are preserved. Values above 1 saturate the
colour scale as in the source; their numerical labels remain unchanged.
The white and hatched cells have no numerical entries in the source.

The drawing style follows EfficiencyTable.plot() in the official systematic
framework: Matplotlib defaults, Seaborn heatmap, mathtext plus/minus,
automatic annotation contrast, original axis labels/title and hatching.
Font sizes, the colourbar gap and outer margins are adjusted for readability.

Requires Python 3, numpy, matplotlib and seaborn.
Run this file to save KaonEff_replotted.png and KaonEff_replotted.pdf beside
the script, or pass --output-stem PATH_WITHOUT_EXTENSION.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle


COS_THETA_EDGES = np.array(
    [-0.866, -0.682, -0.4226, -0.1045, 0.225, 0.5, 0.766, 0.8829, 0.9563]
)
MOMENTUM_EDGES = np.arange(0.0, 4.51, 0.5)  # GeV/c

# Rows run from low to high momentum; columns run from low to high cos(theta).
RATIOS = np.array([
    [1.00, 0.99, 0.99, 1.01, 1.01, 1.01, 1.00, 0.98],
    [0.82, 0.92, 0.99, 0.98, 0.98, 0.98, 0.98, 0.97],
    [0.34, 0.87, 0.98, 0.97, 0.97, 0.97, 0.98, 1.06],
    [0.85, 1.01, 0.98, 0.96, 0.95, 0.97, 1.01, 1.17],
    [0.62, 1.03, 0.97, 0.95, 0.92, 0.94, 0.98, 1.07],
    [0.74, 1.06, 0.94, 0.94, 0.90, 0.92, 0.97, 1.05],
    [0.52, 1.00, 0.92, 0.91, 0.88, 0.90, 0.95, 1.05],
    [1.51, 0.94, 0.88, 0.88, 0.92, 0.90, 0.95, 1.05],
    [np.nan, np.nan, 0.85, 0.92, 0.93, 0.90, 0.94, 1.06],
])
UNCERTAINTIES = np.array([
    [0.01, 0.00, 0.00, 0.00, 0.01, 0.01, 0.01, 0.02],
    [0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.04, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.07, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01],
    [0.08, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.12, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.29, 0.02, 0.01, 0.01, 0.01, 0.00, 0.00, 0.00],
    [0.83, 0.05, 0.01, 0.01, 0.01, 0.00, 0.00, 0.01],
    [np.nan, np.nan, 0.05, 0.03, 0.01, 0.01, 0.01, 0.01],
])

# Preserve the two distinct appearances; do not invent values for these cells.
HATCHED_CELLS = ((8, 1),)
BLANK_CELLS = ((8, 0),)

# The framework uses a 9 x 5 inch canvas. At 0.95\textwidth in
# the thesis, 12.25 pt here corresponds to approximately 8.1 pt.
FIGURE_SIZE_INCHES = (9, 5)
CELL_FONT_SIZE = 12.25
TICK_FONT_SIZE = 12.25
AXIS_FONT_SIZE = 14


def make_figure():
    """Use the framework's two-dimensional heatmap styling."""
    with plt.rc_context():
        # The framework leaves font selection to Matplotlib's defaults.
        # Reset them here so a user's external style does not select Arial.
        plt.rcdefaults()
        plt.rcParams.update({
            "font.size": CELL_FONT_SIZE,
            "xtick.labelsize": TICK_FONT_SIZE,
            "ytick.labelsize": TICK_FONT_SIZE,
            "axes.titlesize": 14,
            "pdf.fonttype": 42,
        })
        fig, ax = plt.subplots(1, 1, figsize=FIGURE_SIZE_INCHES, dpi=120)
        annot_labels = np.asarray([
            f"{value:.2f}" + r"$\pm$" + f"{error:.2f}"
            for value, error in zip(RATIOS.flatten(), UNCERTAINTIES.flatten())
        ]).reshape(RATIOS.shape)
        plot_args = {
            "linewidths": 0.5,
            "fmt": "",
            "cbar_kws": {"fraction": 0.05, "pad": 0.01},
            "vmin": 0,
            # Match the existing thesis plot, whose scale ends at 1.0.
            "vmax": 1,
            "cmap": "plasma",
        }
        sns.heatmap(
            RATIOS, annot=annot_labels, annot_kws={"fontsize": CELL_FONT_SIZE},
            ax=ax, **plot_args,
        )
        ax.set_xticks(np.arange(len(COS_THETA_EDGES)))
        ax.set_xticklabels(np.round(COS_THETA_EDGES, 4))
        ax.set_yticks(np.arange(len(MOMENTUM_EDGES)))
        ax.set_yticklabels(np.round(MOMENTUM_EDGES, 4), rotation=0)
        ax.set_xlabel("cosTheta bins", fontsize=AXIS_FONT_SIZE)
        ax.set_ylabel("p bins", fontsize=AXIS_FONT_SIZE)
        for row, col in HATCHED_CELLS:
            ax.add_patch(Rectangle(
                (col, row), 1, 1, fill=False,
                hatch="////", edgecolor="blue", lw=1,
            ))
        ax.set_title('K ratio table  for cut "kaonID > 0.6"')
        fig.tight_layout(pad=0.5)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-stem", type=Path,
        default=Path(__file__).resolve().with_name("KaonEff_replotted"),
        help="Output path without an extension; its directory must exist.",
    )
    args = parser.parse_args()
    if not args.output_stem.parent.is_dir():
        parser.error(f"Output directory does not exist: {args.output_stem.parent}")

    fig = make_figure()
    for extension in ("png", "pdf"):
        output = Path(f"{args.output_stem}.{extension}")
        with plt.rc_context({"pdf.fonttype": 42, "savefig.bbox": None}):
            fig.savefig(output, dpi=600, facecolor="white", bbox_inches=None)
        print(f"Saved: {output.resolve()}")
    plt.close(fig)
    print("Source: rounded labels transcribed from thesis Fig. K.1.")
    print("70 numerical entries, one blue-hatched cell, one blank cell.")


if __name__ == "__main__":
    main()
