"""Common figure style (close to the MATLAB style of the original figures)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GREEN = (0, 150 / 255, 0)
BLUE = (0, 0, 150 / 255)
RED = (0.85, 0, 0)
MAGENTA = (0.75, 0, 0.75)
CYAN = (0, 0.6, 0.75)
ORANGE = (0.9, 0.5, 0.0)
GREY = (0.4, 0.4, 0.4)
PURPLE = (0.45, 0.2, 0.65)

JOINT_NAMES = ["Base", "Shoulder", "Elbow", "Wrist"]


def setup():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 11,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "legend.fontsize": 9.5,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.linewidth": 1.0,
        "lines.linewidth": 1.5,
        "axes.grid": True,
        "grid.alpha": 0.35,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    })


def save(fig, path):
    fig.savefig(path)
    if path.endswith(".pdf"):
        fig.savefig(path[:-4] + ".png", dpi=150)
    plt.close(fig)
