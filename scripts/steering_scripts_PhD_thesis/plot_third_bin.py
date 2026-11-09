import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EVENT_COUNTS = np.arange(8)
FREQUENCIES = np.array([1, 8, 63, 215, 362, 211, 35, 5])

FIGURE_SIZE_INCHES = (3.14, 3.05)


def make_figure():
    """Return an unnormalized step histogram with publication-size labels."""
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
    }
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE_INCHES)
        fig.subplots_adjust(left=0.23, right=0.97, bottom=0.23, top=0.97)
        ax.stairs(
            FREQUENCIES,
            np.arange(-0.5, 8.0, 1.0),
            baseline=None,
            color="#000099",
            linewidth=1.3,
        )
        ax.set_xlim(-0.5, 7.5)
        ax.set_ylim(0, 380)
        ax.set_xticks(EVENT_COUNTS)
        ax.set_yticks(np.arange(0, 351, 50))
        ax.set_xlabel("Number of data events in the third bin", labelpad=5)
        ax.set_ylabel("Number of retrainings", labelpad=5)
        ax.tick_params(direction="in", length=3, width=0.8, pad=4)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-stem",
        type=Path,
        default=Path(__file__).resolve().with_name("bin_3_replotted"),
        help="Output path without an extension; its directory must exist.",
    )
    args = parser.parse_args()
    if not args.output_stem.parent.is_dir():
        parser.error(f"Output directory does not exist: {args.output_stem.parent}")

    total = int(FREQUENCIES.sum())
    mean = float(np.average(EVENT_COUNTS, weights=FREQUENCIES))
    stddev = float(
        np.sqrt(np.average((EVENT_COUNTS - mean) ** 2, weights=FREQUENCIES))
    )
    print("Source: frequencies reconstructed from the existing PNG.")
    print(f"Entries: {total}; mean: {mean:.6f}; standard deviation: {stddev:.6f}")
    print(f"Exactly one event: {FREQUENCIES[1]}; at most one: {FREQUENCIES[:2].sum()}")

    fig = make_figure()
    for extension in ("png", "pdf"):
        output = Path(f"{args.output_stem}.{extension}")
        # 600 dpi PNG for the existing LaTeX workflow; PDF for vector output.
        with plt.rc_context({"pdf.fonttype": 42, "savefig.bbox": None}):
            fig.savefig(output, dpi=600, facecolor="white", bbox_inches=None)
        print(f"Saved: {output.resolve()}")
    plt.close(fig)


if __name__ == "__main__":
    main()
