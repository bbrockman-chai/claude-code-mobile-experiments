"""Shared chart styling and label layout for the power-analysis figures."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- palette (validated default from the dataviz skill, light mode) ---------
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"

DURATIONS = np.arange(0.5, 5.01, 0.25)
PLAUSIBLE = (0.10, 0.15)  # plausible all-cause U5 mortality effect of CDS/IMCI-type interventions

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "legend.frameon": False,
})


def style_axes(ax, ylabel="Minimum detectable reduction in under-5 mortality (80% power)", durations=DURATIONS):
    ax.set_xlabel("Duration of intervention follow-up (years)")
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.set_xlim(durations[0], durations[-1] + 1.35)
    ax.set_xticks([0.5, 1, 2, 3, 4, 5])
    ax.grid(axis="x", visible=False)


def plausible_band(ax):
    ax.axhspan(PLAUSIBLE[0], PLAUSIBLE[1], color=INK, alpha=0.06, lw=0)
    ax.text(ax.get_xlim()[0] + 0.05, PLAUSIBLE[0] + 0.003,
            "Plausible true effect of CDS on all-cause U5 mortality: 10-15% (IMCI trials)",
            fontsize=8, color=INK2, va="bottom")


_LABELS: dict = {}


def direct_label(ax, x, y, text, color):
    """Queue an end-of-line label; place_labels() lays them out without overlap."""
    _LABELS.setdefault(id(ax), []).append((x, y, text, color))
    ax.plot([x], [y], "o", ms=5, color=color, mec=SURFACE, mew=1.2)


def place_labels(ax, min_gap_frac=0.04):
    """Push queued labels apart vertically so they never overlap; the dot still
    marks the true end point and a hairline leader joins dot to label."""
    items = sorted(_LABELS.pop(id(ax), []), key=lambda t: t[1])
    if not items:
        return
    lo, hi = ax.get_ylim()
    gap = (hi - lo) * min_gap_frac
    xlo, xhi = ax.get_xlim()
    dx = (xhi - xlo) * 0.012
    ys = [it[1] for it in items]
    for i in range(1, len(ys)):                 # push up
        ys[i] = max(ys[i], ys[i - 1] + gap)
    for i in range(len(ys) - 2, -1, -1):        # push down if we overshot the top
        ys[i] = min(ys[i], ys[i + 1] - gap)
    for (x, y, text, color), yl in zip(items, ys):
        if abs(yl - y) > 1e-9:
            ax.plot([x, x + dx], [y, yl], lw=0.8, color=AXIS)
        ax.annotate(text, (x + dx, yl), xytext=(4, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=8.5, color=INK2)


