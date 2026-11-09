"""Draw all four hadron-composition figures in thesis Appendix L.

Weighted yields are taken from the supplied GetFakeFromLepton.c output.
Fractions are normalized within the four displayed categories, as stated
in the thesis captions. The signal-kaon category sum is 22.167345 whereas
the printed all-kaon total is 22.167346; retain the supplied yields without
adding an artificial category for this 0.000001 rounding difference.

The original colours, exploded wedges, start angle, and one-decimal
percentage format are retained. Minor-category labels and percentages are
placed together outside the pie to remain legible without overlapping.
A displayed 0.0% is rounded, not necessarily an exactly zero contribution.

At 0.7--0.75 times the thesis text width, the 10 pt text in this 5-inch
canvas appears at approximately 8.8--9.4 pt. Save with a fixed canvas size.

Run without arguments to save PNG (600 dpi) and vector PDF versions of
signal_Kaon_fake (L.1), background_Kaon_fake (L.2), signal_Pion_fake (L.3),
and background_Pion_fake (L.4) beside this script. Use --output-dir to
select another existing directory, or --chart to draw just one figure.

Requires Python 3, matplotlib, and numpy.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CHARTS = {
    "signal_Kaon_fake": {
        "figure": "L.1",
        "labels": [r"$K$ from $K$", r"$K$ from $\pi$",
                   r"$K$ from $e$", r"$K$ from $\mu$"],
        "reported_total": 22.167346,
        "values": [21.816623, 0.345687, 0.000057, 0.004978],
    },
    "background_Kaon_fake": {
        "figure": "L.2",
        "labels": [r"$K$ from $K$", r"$K$ from $\pi$",
                   r"$K$ from $e$", r"$K$ from $\mu$"],
        "reported_total": 2834.243447,
        "values": [2505.574586, 263.379000, 44.958192, 20.331669],
    },
    "signal_Pion_fake": {
        "figure": "L.3",
        "labels": [r"$\pi$ from $\pi$", r"$\pi$ from $K$",
                   r"$\pi$ from $e$", r"$\pi$ from $\mu$"],
        "reported_total": 9.405912,
        "values": [9.369550, 0.023753, 0.001337, 0.011272],
    },
    "background_Pion_fake": {
        "figure": "L.4",
        "labels": [r"$\pi$ from $\pi$", r"$\pi$ from $K$",
                   r"$\pi$ from $e$", r"$\pi$ from $\mu$"],
        "reported_total": 2137.152490,
        "values": [1982.052356, 22.031521, 26.879752, 106.188861],
    },
}

COLORS = ["lightblue", "deepskyblue", "turquoise", "cornsilk"]
EXPLODE = [0.0, 0.2, 0.35, 0.48]
START_ANGLE = 30
RADIUS = 0.8
FIGURE_SIZE_INCHES = (5.0, 3.4)
FONT_SIZE = 10


def make_figure(name):
    """Draw a pie with readable callouts for the three small categories."""
    chart = CHARTS[name]
    values = np.asarray(chart["values"], dtype=float)
    percentages = 100.0 * values / values.sum()
    with plt.rc_context():
        plt.rcdefaults()
        plt.rcParams.update({"font.size": FONT_SIZE, "pdf.fonttype": 42})
        fig, ax = plt.subplots(figsize=FIGURE_SIZE_INCHES)
        fig.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.97)
        wedges, _ = ax.pie(
            values, explode=EXPLODE, colors=COLORS,
            radius=RADIUS, startangle=START_ANGLE, normalize=True,
        )

        # Preserve the dominant percentage inside the large wedge.
        major_angle = np.deg2rad((wedges[0].theta1 + wedges[0].theta2) / 2.0)
        ax.text(
            0.6 * RADIUS * np.cos(major_angle),
            0.6 * RADIUS * np.sin(major_angle),
            f"{percentages[0]:.1f}%", ha="center", va="center",
        )
        ax.text(0.0, -0.96, chart["labels"][0], ha="center", va="center")

        # Keep minor labels in angular order, with enough vertical spacing.
        # Link each label to its actual (possibly tiny) exploded wedge.
        anchors = []
        for index in range(1, 4):
            wedge = wedges[index]
            angle = np.deg2rad((wedge.theta1 + wedge.theta2) / 2.0)
            x = wedge.center[0] + RADIUS * np.cos(angle)
            y = wedge.center[1] + RADIUS * np.sin(angle)
            anchors.append((y, x, index))
        for rank, (label_y, (wedge_y, wedge_x, index)) in enumerate(zip(
            (0.18, 0.54, 0.90), sorted(anchors)
        )):
            # Bend outside every exploded wedge, so leaders do not cross
            # another category's slice. Stagger the bends to keep them apart.
            elbow_x = 1.40 - 0.035 * rank
            ax.plot(
                [wedge_x, elbow_x, elbow_x, 1.46],
                [wedge_y, wedge_y, label_y, label_y],
                color="0.4", linewidth=0.6,
            )
            ax.text(
                1.50, label_y,
                chart["labels"][index] + f": {percentages[index]:.1f}%",
                ha="left", va="center", fontsize=FONT_SIZE,
            )
        ax.set_xlim(-0.95, 2.65)
        ax.set_ylim(-1.15, 1.15)
        ax.set_aspect("equal")
        ax.set_axis_off()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chart", choices=["all", *CHARTS], default="all")
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent,
        help="Existing output directory; default: this script's directory.",
    )
    args = parser.parse_args()
    if not args.output_dir.is_dir():
        parser.error(f"Output directory does not exist: {args.output_dir}")
    selected = list(CHARTS) if args.chart == "all" else [args.chart]
    for name in selected:
        chart = CHARTS[name]
        values = np.asarray(chart["values"], dtype=float)
        print(
            f"Fig. {chart['figure']}: {name}; "
            f"four-category sum = {values.sum():.6f}; "
            f"reported total = {chart['reported_total']:.6f}"
        )
        print("Percentages:", ", ".join(f"{v:.6f}%" for v in 100 * values / values.sum()))
        fig = make_figure(name)
        for extension in ("png", "pdf"):
            output = args.output_dir / f"{name}.{extension}"
            with plt.rc_context({"pdf.fonttype": 42, "savefig.bbox": None}):
                fig.savefig(output, dpi=600, facecolor="white", bbox_inches=None)
            print(f"Saved: {output.resolve()}")
        plt.close(fig)


if __name__ == "__main__":
    main()
