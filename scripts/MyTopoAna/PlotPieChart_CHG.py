import matplotlib.pyplot as plt
import numpy as np

# [('521=>213|-421|', 11), ('others', 7), ('521=>-423|211|', 6), ('521=>-413|111|211|211|', 6), ('521=>-423|20213|', 6), ('521=>-423|-211|211|211|111|', 4), ('521=>-421|211|', 4), ('521=>-423|213|', 3), ('521=>-413|213|211|', 2), ('521=>-421|113|211|111|', 2), ('521=>-421|113|211|', 2), ('521=>-413|211|211|', 2), ('521=>-411|111|211|211|', 2), ('521=>-421|-211|211|211|', 2), ('521=>-423|211|-211|211|', 2), ('521=>-423|113|211|111|', 1), ('521=>-421|223|211|', 1), ('521=>-411|213|211|', 1), ('521=>20213|-421|', 1), ('521=>-421|113|211|22|', 1), ('521=>-421|431|', 1), ('521=>213|211|-211|211|-211|', 1), ('521=>321|-213|-443|211|', 1), ('521=>-423|-223|211|', 1), ('521=>-423|-221|211|113|', 1), ('521=>-423|321|', 1), ('521=>-423|-211|211|211|111|22|', 1)]

Category = [
r"$B^{+}\rightarrow \bar{D}^{*}(2007)^{0} \mu^{+} \nu_{\mu}$", # 521=>-423|-13|14| 1364
r"$B^{+}\rightarrow \bar{D}^{*}(2007)^{0} e^{+} \nu_{e}$", # 521=>-423|-11|12| 1076
r"$B^{+}\rightarrow \bar{D}^{0} \mu^{+} \nu_{\mu}$",         # 521=>-421|-13|14| 1031
r"$B^{+}\rightarrow \bar{D}^{0} e^{+} \nu_{e}$", # 521=>-421|-11|12| 821
r"$B^{+}\rightarrow \rho(770)^{+} \bar{D}^{0}$", # 521=>213|-421| 624
r"$B^{+}\rightarrow \bar{D}^{*}(2007)^{0} a_{1}(1260)^{+}$",      # 521=>-423|20213| 452
r"$B^{+}\rightarrow \bar{D}^{*}(2007)^{0} \tau^{+} \nu_{\tau}$",   # 521=>-423|-15|16| 403
r"$B^{+}\rightarrow D^{*}(2010)^{-} \pi^{0} \pi^{+} \pi^{+}$",   # 521=>-413|111|211|211| 383
r"$B^{+}\rightarrow \bar{D}^{0} \tau^{+} \nu_{\tau}$",   # 521=>-421|-15|16| 364
r"$B^{+}\rightarrow \bar{D}^{*}(2007)^{0} \pi^{+}$",   # 521=>-423|211| 330
r"$B^{+}\rightarrow \bar{D}^{0} \pi^{+}$",                      # 521=>-421|211| 330
r"$B^{+}\rightarrow \bar{D}^{*}(2007)^{0} \rho(770)^{+}$",   # 521=>-423|213| 317
r"$B^{+}\rightarrow \bar{D}^{*}(2007)^{0} \pi^{-} \pi^{+} \pi^{+} \pi^{0}$",   # 521=>-423|-211|211|211|111| 261
r"$B^{+}\rightarrow \bar{D}^{0} \pi^{-} \pi^{+} \pi^{+}$",   # 521=>-421|-211|211|211| 241
r"$B^{+}\rightarrow \bar{D}^{0} \rho(770)^{0} \pi^{+}$",   # 521=>-421|113|211| 210
"others"]

TotalN = 12812
values = [1364, 1076, 1031, 821, 624, 452, 403, 383, 364, 330, 330, 317, 261, 241, 210]
values.append(TotalN - sum(values))

PANEL_TITLE = '(a) Charged $B\\bar{B}$ background'
OUTPUT_STEM = "CHG_BKG_PIE"

# Draw at the final thesis width (0.98 * 15.92 cm).
FIGURE_WIDTH_CM = 15.6016
FIGURE_HEIGHT_CM = 9.5
FONT_SIZE = 8.2
START_ANGLE = 110


def _spread_labels(positions, minimum_gap):
    """Keep labels in angular order while enforcing a vertical gap."""
    shifted = np.asarray(positions) - np.arange(len(positions)) * minimum_gap
    blocks = []
    for index, position in enumerate(shifted):
        blocks.append([index, index, float(position), 1])
        while len(blocks) > 1 and blocks[-2][2] > blocks[-1][2]:
            right = blocks.pop()
            left = blocks.pop()
            count = left[3] + right[3]
            mean = (left[2] * left[3] + right[2] * right[3]) / count
            blocks.append([left[0], right[1], mean, count])
    result = np.empty(len(positions))
    for first, last, mean, count in blocks:
        result[first:last + 1] = mean
    return result + np.arange(len(positions)) * minimum_gap


def make_figure():
    """Connect every named decay mode directly to its slice."""
    percentages = 100.0 * np.asarray(values, dtype=float) / sum(values)
    width = FIGURE_WIDTH_CM * 72.0 / 2.54
    height = FIGURE_HEIGHT_CM * 72.0 / 2.54
    center_x, center_y = width / 2.0, height / 2.0 - 7.0
    radius = 75.0

    with plt.rc_context({
        "font.family": "DejaVu Sans",
        "font.size": FONT_SIZE,
        "mathtext.fontset": "stix",
        "pdf.fonttype": 42,
    }):
        fig = plt.figure(figsize=(width / 72.0, height / 72.0), dpi=300)
        ax = fig.add_axes([0, 0, 1, 1])
        colors = plt.get_cmap("Spectral")(np.linspace(0, 1, len(values)))
        wedges, _ = ax.pie(
            values,
            colors=colors,
            startangle=START_ANGLE,
            normalize=True,
            center=(center_x, center_y),
            radius=radius,
            wedgeprops={"edgecolor": "white", "linewidth": 0.35},
        )
        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.text(9, height - 10, PANEL_TITLE, ha="left", va="top", fontsize=10)

        callouts = []
        for index, wedge in enumerate(wedges[:-1]):
            angle = np.deg2rad((wedge.theta1 + wedge.theta2) / 2.0)
            cosine, sine = np.cos(angle), np.sin(angle)
            callouts.append({
                "index": index,
                "side": -1 if cosine < 0 else 1,
                "cosine": cosine,
                "sine": sine,
                "y": center_y + (radius + 5) * sine,
            })

        for side in (-1, 1):
            group = sorted(
                (item for item in callouts if item["side"] == side),
                key=lambda item: item["y"],
            )
            if not group:
                continue
            positions = _spread_labels([item["y"] for item in group], 16.0)
            if positions[0] < 25:
                positions += 25 - positions[0]
            if positions[-1] > height - 36:
                positions -= positions[-1] - (height - 36)

            for item, label_y in zip(group, positions):
                index = item["index"]
                cosine, sine = item["cosine"], item["sine"]
                boundary = center_x + side * (radius + 10)
                end_x = 122 if side < 0 else width - 122
                # Exit radially, then route outside the pie in label order.
                points = [
                    (center_x + radius * cosine, center_y + radius * sine),
                    (center_x + (radius + 5) * cosine,
                     center_y + (radius + 5) * sine),
                    (boundary, center_y + (radius + 5) * sine),
                    (end_x, label_y),
                ]
                line_x, line_y = zip(*points)
                ax.plot(line_x, line_y, color="0.3", linewidth=0.55,
                        solid_capstyle="round")
                ax.plot(line_x[0], line_y[0], marker="o",
                        markersize=1.5, color="0.3")
                ax.text(
                    end_x - 3 if side < 0 else end_x + 3,
                    label_y,
                    Category[index] + f"  {percentages[index]:.1f}%",
                    ha="right" if side < 0 else "left",
                    va="center",
                    fontsize=FONT_SIZE,
                )

        # The last entry is the inclusive 'others' category.
        other = wedges[-1]
        angle = np.deg2rad((other.theta1 + other.theta2) / 2.0)
        ax.text(
            center_x + 0.53 * radius * np.cos(angle),
            center_y + 0.53 * radius * np.sin(angle),
            f"Others\n{percentages[-1]:.1f}%",
            ha="center", va="center", color="white",
            fontsize=FONT_SIZE, linespacing=1.3,
        )
    return fig


if __name__ == "__main__":
    with plt.rc_context({"pdf.fonttype": 42}):
        fig = make_figure()
        # Keep the canvas size so the thesis text sizes remain predictable.
        fig.savefig(OUTPUT_STEM + ".pdf", bbox_inches=None)
        fig.savefig(OUTPUT_STEM + ".png", dpi=300, bbox_inches=None)
        plt.close(fig)
