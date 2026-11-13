"""Draw the Figure 5.2 bin-index map as a vector PDF.

Requires Matplotlib; no ROOT files or external LaTeX installation are needed.
By default, bin_index.pdf is saved beside this script. Use --output to choose
another destination. Coordinates follow the physical bin boundaries, and bin
indices increase upwards within each reconstructed X_s mass region.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MASS_EDGES = (0.0, 0.6, 1.0, 2.0)
BDT_EDGES = (
    (0.86, 0.9195, 0.979, 0.986, 0.993, 1.0),  # MR1: indices 1--5
    (0.86, 0.8925, 0.925, 0.95, 0.975, 1.0),  # MR2: indices 6--10
    (0.86, 0.9, 0.94, 0.96, 0.98, 1.0),      # MR3: indices 11--15
)


def draw_bin_index(output: Path) -> None:
    """Preserve the unequal bin widths and heights when drawing the map."""
    style = {
        "font.family": "STIXGeneral",
        "mathtext.fontset": "stix",
        "text.usetex": False,
        "pdf.fonttype": 42,
        "axes.linewidth": 1.1,
    }
    with plt.rc_context(style):
        # Extra height leaves space around the labels in the narrow bins 3--5.
        fig, ax = plt.subplots(figsize=(6.4, 5.8))
        fig.subplots_adjust(left=0.16, right=0.97, bottom=0.18, top=0.97)

        ax.set_xlim(MASS_EDGES[0], MASS_EDGES[-1])
        ax.set_ylim(BDT_EDGES[0][0], BDT_EDGES[0][-1])
        ax.vlines(MASS_EDGES[1:-1], 0.86, 1.0, colors="black", linewidth=1.1)

        for region, (mass_low, mass_high) in enumerate(
            zip(MASS_EDGES[:-1], MASS_EDGES[1:])
        ):
            edges = BDT_EDGES[region]
            ax.hlines(edges[1:-1], mass_low, mass_high,
                      colors="black", linewidth=1.1)
            for local_bin, (bdt_low, bdt_high) in enumerate(
                zip(edges[:-1], edges[1:]), start=1
            ):
                ax.text(
                    (mass_low + mass_high) / 2,
                    (bdt_low + bdt_high) / 2,
                    str(5 * region + local_bin),
                    ha="center", va="center", fontsize=14,
                )

        ax.set_xticks(MASS_EDGES, labels=("0.0", "0.6", "1.0", "2.0"))
        ax.set_yticks((0.86, 1.0), labels=("0.86", "1.0"))
        ax.tick_params(axis="both", which="both", length=0,
                       pad=6, labelsize=14)
        ax.set_xlabel(r"$M_{X_s}^{\mathrm{reco}}$ [GeV/$c^2$]",
                      fontsize=18, labelpad=10)
        ax.set_ylabel(r"$\mathcal{O}_{\mathrm{BDT}}$",
                      fontsize=18, labelpad=10)

        output = output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, format="pdf", bbox_inches="tight",
                    pad_inches=0.06, facecolor="white")
        plt.close(fig)
    print(f"Saved {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().with_name("bin_index.pdf"),
        help="Output PDF path (default: bin_index.pdf beside this script).",
    )
    args = parser.parse_args()
    draw_bin_index(args.output)


if __name__ == "__main__":
    main()
