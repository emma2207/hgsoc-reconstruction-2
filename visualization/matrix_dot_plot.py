
import numpy as np
import matplotlib.pyplot as plt
from parameters import TOOLS


def matrix_dot_plot():

    criteria = [
        "Pseudobulk performance", "HGSOC performance",
        "Computational cost", "Usability"
    ]

    # 1 = black, 0.5 = half-filled, 0 = white
    # -1 = gray
    data = np.array([
        [1, 1, 0, .5],  # Conpair
        [0, .5, .5, .5],  # Crosscheck
        [.5, 0, .5, .5],  # HYSYS
        [1, 1, .5, 1],  # NGSCheckMate
        [1, .5, 1, 1],  # ntsm
        [.5, .5, .5, .5],  # OmicsPrint
        [1, 1, .5, .5],  # Peddy
        [0, 0, 1, 1],  # Somalier
        [.5, 0, .5, .5],  # TimeAttack
        [1, .5, .5, .5],  # Vireo
    ])

    fig, ax = plt.subplots(figsize=(7, 7))

    # Horizontal guide lines
    for y in range(len(TOOLS)):
        ax.axhline(y, color="#dddddd", lw=1.5, zorder=0)

    # Draw circles
    for i in range(len(TOOLS)):
        for j in range(len(criteria)):
            value = data[i, j]

            if value == 1:
                fillstyle = "full"
                color = "black"
            elif value == 0:
                fillstyle = "full"
                color = "white"
            elif value == 0.5:
                fillstyle = "bottom"
                color = "black"
            else:
                fillstyle = "full"
                color = "lightgray"

            ax.plot(
                j, i,
                marker="o",
                markersize=26,
                linestyle="None",
                markerfacecolor=color,
                markerfacecoloralt="white",
                markeredgecolor=(
                    "lightgray" if value == -1 else "black"
                ),
                markeredgewidth=1.5,
                fillstyle=fillstyle,
                zorder=2,
            )

    # Axis formatting
    ax.set_xticks(range(len(criteria)), criteria)
    ax.set_yticks(range(len(TOOLS)), TOOLS)

    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=90, fontsize=13)
    plt.setp(ax.get_yticklabels(), fontsize=13)

    ax.tick_params(axis="both", length=0, pad=12)
    ax.set_xlim(-0.6, len(criteria) - 0.5)
    ax.set_ylim(len(TOOLS) - 0.5, -0.5)

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.savefig("../benchmarking/data/figures/matrix_dot_plot.png", dpi=300)
    plt.savefig("../benchmarking/data/figures/matrix_dot_plot.pdf", dpi=300)


if __name__ == "__main__":
    matrix_dot_plot()
