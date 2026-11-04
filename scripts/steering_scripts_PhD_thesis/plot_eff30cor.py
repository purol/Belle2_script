#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Draw a data/MC ratio heatmap from eff30cor.csv."
    )
    parser.add_argument(
        "input_csv",
        nargs="?",
        default="eff30cor.csv",
        help="Input CSV file (default: eff30cor.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="eff30cor.pdf",
        help="Output image file (default: eff30cor.pdf)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_csv)
    output_path = Path(args.output)

    df = pd.read_csv(input_path)

    required_columns = {
        "p_min",
        "p_max",
        "cosTheta_min",
        "cosTheta_max",
        "data_MC_ratio",
        "data_MC_ratio_err",
    }
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(
            "Missing required columns: " + ", ".join(sorted(missing_columns))
        )

    # Obtain the bin boundaries directly from the CSV file.
    p_edges = np.array(
        sorted(set(df["p_min"].tolist()) | set(df["p_max"].tolist())),
        dtype=float,
    )
    cos_theta_edges = np.array(
        sorted(
            set(df["cosTheta_min"].tolist())
            | set(df["cosTheta_max"].tolist())
        ),
        dtype=float,
    )

    n_p_bins = len(p_edges) - 1
    n_cos_theta_bins = len(cos_theta_edges) - 1

    ratio_array = np.full((n_p_bins, n_cos_theta_bins), np.nan)
    annotation_array = np.full(
        (n_p_bins, n_cos_theta_bins), "", dtype=object
    )

    for row in df.itertuples(index=False):
        p_matches = np.where(np.isclose(p_edges[:-1], row.p_min))[0]
        cos_theta_matches = np.where(
            np.isclose(cos_theta_edges[:-1], row.cosTheta_min)
        )[0]

        if len(p_matches) != 1 or len(cos_theta_matches) != 1:
            raise ValueError(
                "Could not identify a unique bin for "
                f"p_min={row.p_min}, cosTheta_min={row.cosTheta_min}"
            )

        p_index = p_matches[0]
        cos_theta_index = cos_theta_matches[0]

        ratio_array[p_index, cos_theta_index] = row.data_MC_ratio
        annotation_array[p_index, cos_theta_index] = (
            f"{row.data_MC_ratio:.2f}"
            f"$\\pm${row.data_MC_ratio_err:.2f}"
        )

    fig, ax = plt.subplots(figsize=(11.0, 6.0), dpi=120)

    sns.heatmap(
        ratio_array,
        annot=annotation_array,
        fmt="",
        linewidths=0.5,
        vmin=0.0,
        vmax=1.0,
        cmap="plasma",
        cbar_kws={"fraction": 0.05, "label": ""},
        ax=ax,
        annot_kws={"fontsize": 15}
    )

    # Heatmap cells are centered at 0.5, 1.5, ...; boundary labels are placed
    # at integer coordinates to reproduce the bin-edge style of the reference.
    ax.set_xticks(np.arange(len(cos_theta_edges)))
    ax.set_xticklabels([f"{value:g}" for value in cos_theta_edges])
    ax.set_yticks(np.arange(len(p_edges)))
    ax.set_yticklabels([f"{value:g}" for value in p_edges], rotation=0)

    ax.tick_params(axis="both", labelsize=15)
    ax.collections[0].colorbar.ax.tick_params(labelsize=15)

    ax.set_xlabel(r"$\cos\theta$ bins", fontsize=15)
    ax.set_ylabel(r"$p$ bins", fontsize=15)
    ax.set_title(r"$\pi^{0}$ Efficiency Data/MC Ratio", fontsize=15)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved plot to {output_path}")


if __name__ == "__main__":
    main()
