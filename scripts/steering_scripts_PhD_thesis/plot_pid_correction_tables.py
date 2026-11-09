"""Replot the remaining PID tables (thesis Figs. K.2--K.12).

All numerical labels were transcribed from the PNG figures under
phd_thesis/picture/systematic, NOT from the original calibration CSV files.
Their displayed precision, signs, NaN uncertainty, empty cells, and blue
hatching are preserved. Do not use these rounded values for PID correction.

A row token is value/error for symmetric errors, or value/up/down for
asymmetric errors. B denotes a blank cell; H denotes a blue-hatched cell.
The asymmetric uncertainties retain the source superscript/subscript style,
and all muon tables retain eight columns, as requested.

The drawing style follows the systematic framework and existing PID plots:
Matplotlib defaults (DejaVu Sans), sns.heatmap, mathtext, automatic text
contrast, plasma colours, and blue hatching. All colour scales remain 0--1,
as in the source PNGs; the numerical labels retain ratios above 1.

Run without options to draw all eleven tables. Use --table NAME for one
table, and --output-dir DIRECTORY to select an existing output directory.
Outputs use *_replotted.png (600 dpi) and *_replotted.pdf; source figures
and the existing plot_kaon_eff_table.py are not touched.

Recommended thesis widths: 0.95\\textwidth for all tables except Pionmiss
(Fig. K.8), which needs 1.0\\textwidth for its longer numerical labels.
Muon tables keep eight columns and smaller superscript/subscript errors.

Requires Python 3, numpy, matplotlib, and seaborn.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle


# The rows follow increasing momentum; columns follow increasing angle.
TABLES = {
    "Kaonmis": {
        "figure": "K.2",
        "kind": "symmetric",
        "title": "pi ratio table  for cut \"kaonID > 0.6\"",
        "source": "Kaonmis.png",
        "rows": [
            "1.97/2.77 1.51/0.11 1.59/0.09 2.16/0.24 1.42/0.27 1.25/0.09 1.27/0.28 2.08/1.08",
            "0.46/0.06 0.92/0.03 1.41/0.04 1.48/0.03 1.46/0.03 1.34/0.03 0.92/0.03 0.52/0.04",
            "0.25/0.07 0.54/0.03 1.01/0.04 1.12/0.03 1.33/0.03 1.27/0.02 0.78/0.03 0.58/0.03",
            "0.37/0.14 1.57/0.03 1.13/0.04 1.14/0.03 1.32/0.02 1.50/0.02 1.27/0.03 0.97/0.05",
            "0.29/0.24 3.11/0.06 1.32/0.05 1.20/0.04 1.18/0.02 1.43/0.02 1.71/0.04 1.23/0.08",
            "0.44/0.37 2.83/0.04 1.19/0.05 1.16/0.04 0.96/0.02 1.31/0.02 1.49/0.03 1.05/0.09",
            "0.76/1.14 1.98/0.05 1.04/0.06 0.91/0.04 0.84/0.03 1.08/0.02 1.07/0.03 1.53/0.11",
            "H 1.71/0.35 1.03/0.14 0.87/0.06 0.77/0.04 0.96/0.03 0.88/0.04 1.45/0.13",
            "B H 1.57/2.31 0.93/0.15 0.76/0.07 1.07/0.05 0.67/0.05 1.42/0.15"
        ],
        "x_edges": [
            -0.866,
            -0.682,
            -0.4226,
            -0.1045,
            0.225,
            0.5,
            0.766,
            0.8829,
            0.9563
        ],
        "y_edges": [
            0,
            0.5,
            1,
            1.5,
            2,
            2.5,
            3,
            3.5,
            4,
            4.5
        ]
    },
    "E_plus_kaonID": {
        "figure": "K.3",
        "kind": "asymmetric",
        "title": "e$^{+}$ ratio table for cut \"kaonID>0.6\"",
        "source": "E_plus_kaonID.png",
        "rows": [
            "0.62/0.05/0.05 1.20/0.04/0.06 1.00/0.16/0.32 0.91/0.20/0.45 1.19/0.09/0.20 0.86/0.04/0.06",
            "0.66/0.01/0.01 1.11/0.03/0.05 0.56/0.35/0.44 0.20/0.48/0.25 1.40/0.06/0.14 1.12/0.01/0.03",
            "0.46/0.09/0.03 B 0.58/0.28/0.22 B 1.78/0.15/0.09 0.79/0.04/0.03",
            "B B 1.30/1.23/0.68 B 0.00/0.00/0.00 B",
            "0.59/0.44/0.62 0.72/0.12/0.10 B B B B",
            "1.82/2.95/1.59 B B B B B"
        ],
        "x_edges": [
            0.22,
            0.56,
            1.13,
            1.57,
            1.88,
            2.23,
            2.71
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5,
            3
        ]
    },
    "E_minus_kaonID": {
        "figure": "K.4",
        "kind": "asymmetric",
        "title": "e$^{-}$ ratio table for cut \"kaonID>0.6\"",
        "source": "E_minus_kaonID.png",
        "rows": [
            "0.83/0.07/0.07 B B B 0.80/0.24/0.61 0.92/0.04/0.08",
            "0.48/0.08/0.14 0.70/0.20/0.26 B B 1.14/0.09/0.19 0.95/0.02/0.08",
            "B B 1.04/0.04/0.06 1.06/0.09/0.20 1.54/0.05/0.08 0.84/0.02/0.02",
            "B 1.35/0.12/0.15 1.68/0.56/0.54 0.15/0.78/0.30 0.00/0.00/0.00 B",
            "1.91/0.15/0.15 0.81/0.08/0.08 B B B B",
            "0.00/0.00/0.00 B B B B B"
        ],
        "x_edges": [
            0.22,
            0.56,
            1.13,
            1.57,
            1.88,
            2.23,
            2.71
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5,
            3
        ]
    },
    "MU_plus_kaonID": {
        "figure": "K.5",
        "kind": "asymmetric_muon",
        "title": "$\\mu^{+}$ ratio table for cut \"kaonID>0.6\"",
        "source": "MU_plus_kaonID.png",
        "rows": [
            "0.63/0.10/0.07 2.28/0.15/0.20 1.88/0.08/0.11 2.29/0.06/0.09 2.49/0.08/0.13 1.85/0.05/0.08 1.06/0.02/0.03 0.96/0.03/0.05",
            "0.43/0.04/0.04 3.47/0.23/0.32 3.35/0.23/0.37 3.45/0.21/0.34 4.08/0.21/0.38 3.22/0.15/0.23 1.34/0.08/0.09 0.78/0.07/0.08",
            "0.46/0.11/0.10 2.34/0.10/0.14 2.47/0.11/0.17 2.68/0.17/0.25 2.71/0.28/0.44 1.78/0.13/0.16 1.07/0.04/0.04 0.88/0.05/0.05",
            "0.67/0.06/0.13 3.97/0.60/0.76 3.44/0.55/0.70 15.18/11.13/7.93 2.17/2.29/1.49 0.00/0.00/0.00 2.61/4.18/2.21 0.00/0.00/0.00",
            "1.05/0.30/0.53 3.41/2.66/2.01 0.00/0.00/0.00 B B B B B"
        ],
        "x_edges": [
            0.4,
            0.64,
            0.82,
            1.16,
            1.46,
            1.78,
            2.13,
            2.22,
            2.6
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5
        ]
    },
    "MU_minus_kaonID": {
        "figure": "K.6",
        "kind": "asymmetric_muon",
        "title": "$\\mu^{-}$ ratio table for cut \"kaonID>0.6\"",
        "source": "MU_minus_kaonID.png",
        "rows": [
            "0.74/0.07/0.15 2.13/0.14/0.16 1.79/0.07/0.10 2.27/0.07/0.10 2.38/0.08/0.15 1.85/0.05/0.08 1.21/0.03/0.04 1.03/0.03/0.07",
            "0.81/0.06/0.06 3.06/0.18/0.20 3.31/0.20/0.24 3.63/0.20/0.30 3.97/0.20/0.36 3.77/0.17/0.23 1.74/0.10/0.11 1.23/0.12/0.12",
            "0.52/0.08/0.27 1.96/0.10/0.12 2.39/0.12/0.16 2.27/0.18/0.29 2.24/0.26/0.49 1.85/0.16/0.25 0.87/0.04/0.06 0.61/0.04/0.08",
            "0.72/0.06/0.09 4.91/0.85/0.82 4.71/0.94/0.95 0.00/0.00/0.00 B 0.00/0.00/0.00 0.00/0.00/0.00 B",
            "2.89/1.20/1.04 3.18/2.06/1.52 0.00/0.00/0.00 B B B B B"
        ],
        "x_edges": [
            0.4,
            0.64,
            0.82,
            1.16,
            1.46,
            1.78,
            2.13,
            2.22,
            2.6
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5
        ]
    },
    "PionEff": {
        "figure": "K.7",
        "kind": "symmetric",
        "title": "pi ratio table  for cut \"pionID > 0.6\"",
        "source": "PionEff.png",
        "rows": [
            "0.85/0.04 0.95/0.01 0.95/0.01 0.94/0.01 0.98/0.01 0.96/0.01 0.88/0.01 0.88/0.02",
            "1.16/0.01 1.11/0.01 0.95/0.00 0.97/0.00 0.99/0.00 1.00/0.00 0.96/0.00 0.96/0.01",
            "0.91/0.03 0.98/0.01 0.95/0.00 0.95/0.00 0.96/0.00 0.97/0.00 0.99/0.00 1.05/0.01",
            "0.38/0.11 0.93/0.01 0.96/0.00 0.95/0.00 0.95/0.00 0.95/0.00 0.96/0.00 1.04/0.01",
            "0.22/0.11 0.77/0.01 0.93/0.00 0.91/0.00 0.92/0.00 0.91/0.00 0.94/0.00 0.99/0.01",
            "0.27/0.12 0.71/0.01 0.94/0.00 0.92/0.00 0.94/0.00 0.89/0.00 0.94/0.00 0.99/0.01",
            "0.19/0.33 0.72/0.03 0.95/0.01 0.94/0.01 0.95/0.01 0.90/0.00 0.96/0.00 0.97/0.01",
            "-0.00/nan 0.77/0.20 0.94/0.02 0.96/0.01 0.97/0.01 0.90/0.01 0.98/0.01 0.97/0.01",
            "B B 0.60/0.66 0.94/0.04 0.99/0.02 0.91/0.01 1.05/0.01 1.00/0.01"
        ],
        "x_edges": [
            -0.866,
            -0.682,
            -0.4226,
            -0.1045,
            0.225,
            0.5,
            0.766,
            0.8829,
            0.9563
        ],
        "y_edges": [
            0,
            0.5,
            1,
            1.5,
            2,
            2.5,
            3,
            3.5,
            4,
            4.5
        ]
    },
    "Pionmiss": {
        "figure": "K.8",
        "kind": "symmetric",
        "title": "K ratio table  for cut \"pionID > 0.6\"",
        "source": "Pionmiss.png",
        "rows": [
            "H 2.16/8.37 1.23/1.12 0.59/16.03 0.81/17.25 H 0.50/3.82 0.32/4.42",
            "2.60/0.12 1.74/0.08 1.00/0.06 1.16/0.04 1.19/0.04 1.35/0.08 1.31/0.20 2.04/1.82",
            "1.66/0.05 1.03/0.03 0.96/0.03 1.03/0.02 1.15/0.02 1.24/0.02 1.38/0.04 1.93/0.17",
            "0.43/0.14 0.71/0.03 0.96/0.03 1.11/0.02 1.21/0.02 1.16/0.02 1.42/0.03 1.71/0.13",
            "0.21/0.18 0.66/0.03 1.21/0.03 1.18/0.02 1.19/0.01 1.25/0.02 1.45/0.03 1.63/0.12",
            "0.22/0.23 0.66/0.03 1.45/0.03 1.17/0.02 1.12/0.01 1.14/0.01 1.52/0.03 1.53/0.14",
            "0.18/0.37 0.66/0.04 1.47/0.03 1.14/0.02 1.09/0.01 1.15/0.02 1.62/0.03 1.57/0.18",
            "0.16/1.02 0.64/0.10 1.66/0.05 1.13/0.03 0.99/0.02 1.15/0.02 1.62/0.03 1.63/0.18",
            "B H 1.59/0.18 0.93/0.06 0.93/0.04 1.17/0.04 1.66/0.04 1.41/0.10"
        ],
        "x_edges": [
            -0.866,
            -0.682,
            -0.4226,
            -0.1045,
            0.225,
            0.5,
            0.766,
            0.8829,
            0.9563
        ],
        "y_edges": [
            0,
            0.5,
            1,
            1.5,
            2,
            2.5,
            3,
            3.5,
            4,
            4.5
        ]
    },
    "E_plus_pionID": {
        "figure": "K.9",
        "kind": "asymmetric",
        "title": "e$^{+}$ ratio table for cut \"pionID>0.6\"",
        "source": "E_plus_pionID.png",
        "rows": [
            "0.88/0.09/0.09 1.28/0.29/0.12 0.64/0.49/0.20 B 2.05/0.44/0.12 0.78/0.09/0.06",
            "2.13/0.36/0.38 0.64/0.79/0.47 B B 2.73/0.86/0.42 0.45/0.23/0.07",
            "1.53/0.70/0.60 B B B 1.48/1.03/0.80 B",
            "0.58/0.40/0.17 B B B 0.34/1.53/1.35 B",
            "0.13/0.17/0.07 0.05/0.16/0.11 B B B B",
            "B B B B B B"
        ],
        "x_edges": [
            0.22,
            0.56,
            1.13,
            1.57,
            1.88,
            2.23,
            2.71
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5,
            3
        ]
    },
    "E_minus_pionID": {
        "figure": "K.10",
        "kind": "asymmetric",
        "title": "e$^{-}$ ratio table for cut \"pionID>0.6\"",
        "source": "E_minus_pionID.png",
        "rows": [
            "0.87/0.07/0.07 1.25/0.10/0.02 0.89/0.20/0.05 0.43/0.51/0.12 1.65/0.23/0.05 0.30/0.02/0.01",
            "1.95/0.27/0.08 1.13/0.35/0.09 0.60/0.68/0.12 0.52/0.73/0.14 2.99/0.37/0.07 0.29/0.13/0.06",
            "1.55/0.55/0.11 0.72/0.14/0.11 1.11/0.93/1.38 1.21/0.35/0.61 3.17/0.13/0.17 B",
            "0.93/0.04/0.06 0.66/0.10/0.10 1.65/1.63/2.30 B 3.44/2.66/1.99 B",
            "0.30/0.02/0.02 0.31/0.02/0.02 B B B B",
            "0.00/0.17/0.14 B B B B B"
        ],
        "x_edges": [
            0.22,
            0.56,
            1.13,
            1.57,
            1.88,
            2.23,
            2.71
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5,
            3
        ]
    },
    "MU_plus_pionID": {
        "figure": "K.11",
        "kind": "asymmetric_muon",
        "title": "$\\mu^{+}$ ratio table for cut \"pionID>0.6\"",
        "source": "MU_plus_pionID.png",
        "rows": [
            "1.54/0.10/0.08 1.67/0.11/0.15 1.63/0.09/0.14 1.85/0.04/0.08 1.59/0.04/0.10 1.81/0.02/0.04 2.31/0.03/0.03 2.09/0.03/0.05",
            "1.02/0.02/0.02 1.24/0.02/0.03 1.37/0.06/0.11 1.57/0.13/0.25 1.41/0.12/0.29 1.84/0.05/0.11 1.82/0.02/0.04 2.20/0.04/0.06",
            "1.76/0.11/0.10 1.87/0.05/0.08 2.37/0.05/0.09 2.33/0.08/0.15 2.13/0.10/0.19 1.96/0.05/0.06 2.52/0.07/0.07 1.61/0.07/0.08",
            "1.56/0.03/0.06 2.86/0.10/0.16 2.38/0.12/0.18 2.09/0.20/0.29 1.83/0.41/0.46 2.94/0.94/0.80 1.80/3.06/1.67 B",
            "1.26/0.05/0.10 2.02/0.22/0.23 1.48/1.52/1.00 B B B B B"
        ],
        "x_edges": [
            0.4,
            0.64,
            0.82,
            1.16,
            1.46,
            1.78,
            2.13,
            2.22,
            2.6
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5
        ]
    },
    "MU_minus_pionID": {
        "figure": "K.12",
        "kind": "asymmetric_muon",
        "title": "$\\mu^{-}$ ratio table for cut \"pionID>0.6\"",
        "source": "MU_minus_pionID.png",
        "rows": [
            "1.60/0.08/0.16 1.92/0.11/0.12 1.79/0.09/0.12 1.90/0.04/0.09 1.58/0.04/0.12 1.82/0.02/0.04 2.35/0.03/0.04 2.16/0.03/0.06",
            "1.04/0.02/0.02 1.30/0.01/0.02 1.45/0.06/0.07 1.62/0.13/0.22 1.41/0.12/0.31 1.87/0.05/0.11 1.82/0.02/0.04 2.03/0.04/0.06",
            "1.82/0.08/0.24 1.99/0.05/0.06 2.39/0.05/0.08 2.27/0.08/0.15 1.89/0.09/0.20 1.81/0.04/0.06 2.58/0.06/0.07 1.96/0.08/0.09",
            "1.54/0.03/0.05 2.63/0.10/0.11 2.36/0.13/0.16 1.92/0.20/0.26 1.69/0.40/0.43 1.17/0.23/0.22 5.48/8.39/4.36 1.21/2.20/1.22",
            "1.24/0.06/0.07 2.31/0.27/0.25 1.17/0.82/0.61 B B B B B"
        ],
        "x_edges": [
            0.4,
            0.64,
            0.82,
            1.16,
            1.46,
            1.78,
            2.13,
            2.22,
            2.6
        ],
        "y_edges": [
            0.4,
            0.5,
            1,
            1.5,
            2,
            2.5
        ]
    }
}


def parse_table(table):
    """Return values, exact displayed annotations, and hatched cell indices."""
    rows = [row.split() for row in table["rows"]]
    expected = (len(table["y_edges"]) - 1, len(table["x_edges"]) - 1)
    if len(rows) != expected[0] or any(len(row) != expected[1] for row in rows):
        raise ValueError(f"Invalid table dimensions for Fig. {table['figure']}")
    values = np.full(expected, np.nan)
    labels = np.full(expected, "", dtype=object)
    hatched = []
    symmetric = table["kind"] == "symmetric"
    for i, row in enumerate(rows):
        for j, token in enumerate(row):
            if token in ("B", "H"):
                if token == "H":
                    hatched.append((i, j))
                continue
            parts = token.split("/")
            if len(parts) != (2 if symmetric else 3):
                raise ValueError(f"Invalid label in Fig. {table['figure']}: {token}")
            values[i, j] = float(parts[0])
            if symmetric:
                labels[i, j] = parts[0] + r"$\pm$" + parts[1]
            else:
                labels[i, j] = parts[0] + "$^{+" + parts[1] + "}_{-" + parts[2] + "}$"
    return values, labels, hatched


def plot_settings(name):
    """Keep original canvas proportions; enlarge text for thesis use."""
    kind = TABLES[name]["kind"]
    # Two errors in K.8 have an extra digit (16.03 and 17.25). Use the
    # full thesis text width for this table so its 11.5 pt labels stay >=8 pt.
    if name == "Pionmiss":
        return (9, 5), 11.5, 12.25, 14
    if kind == "symmetric":
        return (9, 5), 12.25, 12.25, 14
    if kind == "asymmetric_muon":
        return (10, 5), 14.25, 14, 15
    return (10, 5), 16, 14, 15


def make_figure(name):
    """Draw one table with the source heatmap style and one annotation per cell."""
    table = TABLES[name]
    values, labels, hatched = parse_table(table)
    figsize, cell_font, tick_font, axis_font = plot_settings(name)
    with plt.rc_context():
        plt.rcdefaults()
        plt.rcParams.update({
            "font.size": cell_font,
            "xtick.labelsize": tick_font,
            "ytick.labelsize": tick_font,
            "axes.titlesize": axis_font,
            "pdf.fonttype": 42,
        })
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=120)
        sns.heatmap(
            values, annot=labels, annot_kws={"fontsize": cell_font},
            ax=ax, linewidths=0.5, fmt="",
            cbar_kws={"fraction": 0.05, "pad": 0.005 if name == "Pionmiss" else 0.01},
            vmin=0, vmax=1, cmap="plasma",
        )
        ax.set_xticks(np.arange(len(table["x_edges"])))
        ax.set_xticklabels(np.round(table["x_edges"], 4))
        ax.set_yticks(np.arange(len(table["y_edges"])))
        ax.set_yticklabels(np.round(table["y_edges"], 4), rotation=0)
        if table["kind"] == "symmetric":
            ax.set_xlabel("cosTheta bins", fontsize=axis_font)
            ax.set_ylabel("p bins", fontsize=axis_font)
        else:
            ax.set_xlabel(r"$\theta$ [rad]", fontsize=axis_font)
            ax.set_ylabel(r"$p$ [GeV/$c$]", fontsize=axis_font)
        for row, col in hatched:
            ax.add_patch(Rectangle(
                (col, row), 1, 1, fill=False,
                hatch="////", edgecolor="blue", lw=1,
            ))
        ax.set_title(table["title"])
        fig.tight_layout(pad=0.3 if name == "Pionmiss" else 0.5)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--table", choices=["all", *TABLES], default="all",
        help="Table to draw; default: all eleven remaining PID tables.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent,
        help="Existing output directory; default: this script's directory.",
    )
    args = parser.parse_args()
    if not args.output_dir.is_dir():
        parser.error(f"Output directory does not exist: {args.output_dir}")
    selected = list(TABLES) if args.table == "all" else [args.table]
    for name in selected:
        fig = make_figure(name)
        for extension in ("png", "pdf"):
            output = args.output_dir / f"{name}_replotted.{extension}"
            with plt.rc_context({"pdf.fonttype": 42, "savefig.bbox": None}):
                fig.savefig(output, dpi=600, facecolor="white", bbox_inches=None)
            print(f"Saved Fig. {TABLES[name]['figure']}: {output.resolve()}")
        plt.close(fig)
        width = "1.0" if name == "Pionmiss" else "0.95"
        print("Recommended LaTeX width: " + width + r"\textwidth")
    print("Source: rounded numerical labels transcribed from thesis plots.")


if __name__ == "__main__":
    main()
