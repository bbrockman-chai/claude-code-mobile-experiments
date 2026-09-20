"""Generate the figures and tables for the Lesotho CDS under-5 mortality power analysis.

Run:  python3 run_analysis.py
Outputs go to figures/ and results/.
"""
from __future__ import annotations

import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import (Assumptions, mdes_parallel, mdes_stepped_wedge, mdes_did,
                   power_parallel, power_did, power_stepped_wedge, expected_events,
                   residual_kb2)

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
RES = os.path.join(HERE, "results")
os.makedirs(FIG, exist_ok=True)
os.makedirs(RES, exist_ok=True)

# --- palette (validated default from the dataviz skill, light mode) ---------
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"

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

BASE = Assumptions()            # see README for the provenance of every default
DURATIONS = np.arange(0.5, 5.01, 0.25)
PLAUSIBLE = (0.10, 0.15)        # plausible all-cause U5 mortality effect of CDS/IMCI-type interventions


def style_axes(ax, ylabel="Minimum detectable reduction in under-5 mortality (80% power)"):
    ax.set_xlabel("Duration of intervention follow-up (years)")
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.set_xlim(DURATIONS[0], DURATIONS[-1] + 1.35)
    ax.set_xticks([0.5, 1, 2, 3, 4, 5])
    ax.grid(axis="x", visible=False)


def plausible_band(ax):
    ax.axhspan(PLAUSIBLE[0], PLAUSIBLE[1], color=INK, alpha=0.06, lw=0)
    ax.text(DURATIONS[0] + 0.05, PLAUSIBLE[0] + 0.003,
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


# ============================================================================
# Figure 1: randomised designs, MDES vs duration
# ============================================================================
designs_rct = [
    ("A. Parallel CRT, 30 vs 30, no baseline",
     lambda T: mdes_parallel(BASE, T)),
    ("B. Parallel CRT, 30 vs 90 controls, no baseline",
     lambda T: mdes_parallel(BASE.with_(c0=90), T)),
    ("C. Parallel CRT, 30 vs 30, adjusted for 3-yr DHIS2 baseline",
     lambda T: mdes_parallel(BASE.with_(baseline_years=3), T)),
    ("D. Parallel CRT, 30 vs 90, adjusted for 3-yr DHIS2 baseline",
     lambda T: mdes_parallel(BASE.with_(c0=90, baseline_years=3), T)),
    ("E. Stepped wedge, 30 facilities, 5 steps (no external controls)",
     lambda T: mdes_stepped_wedge(BASE, T, steps=5)),
]
colors_rct = [SERIES[0], SERIES[1], SERIES[2], SERIES[3], SERIES[4]]
short_rct = ["A  30 v 30", "B  30 v 90", "C  30 v 30 + baseline", "D  30 v 90 + baseline", "E  stepped wedge"]

table_rows = []
fig, ax = plt.subplots(figsize=(10, 6.2))
plausible_band(ax)
for (label, fn), col, short in zip(designs_rct, colors_rct, short_rct):
    y = np.array([fn(T) for T in DURATIONS])
    ax.plot(DURATIONS, y, lw=2, color=col, label=label)
    direct_label(ax, DURATIONS[-1], y[-1], short, col)
    for T, v in zip(DURATIONS, y):
        table_rows.append((label, T, v))
style_axes(ax)
ax.set_ylim(0, 0.45)
place_labels(ax)
ax.set_title("What effect on under-5 mortality could a 30-facility CDS trial in Lesotho detect?\n"
             "Randomised designs; high-burden catchments (U5MR ~65/1,000, ~250 births/facility/yr)",
             loc="left", fontsize=11.5)
ax.legend(loc="upper right", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig1_randomised_designs_mdes_vs_duration.png"), dpi=170)
plt.close(fig)

# ============================================================================
# Figure 2: adding the non-randomised secondary analyses
# ============================================================================
designs_sec = [
    ("D. Randomised: 30 vs 90, adjusted for 3-yr baseline (best RCT above)",
     lambda T: mdes_parallel(BASE.with_(c0=90, baseline_years=3), T), SERIES[3], "D  best RCT"),
    ("F. Non-randomised DiD vs 150 DHIS2 facilities, 3-yr baseline (DHIS2 outcome, 70% complete)",
     lambda T: mdes_did(BASE.with_(completeness=0.7), T, B=3, n_comparison=150), SERIES[6], "F  DiD, 3-yr baseline"),
    ("G. Non-randomised DiD vs 150 DHIS2 facilities, 6-yr baseline (DHIS2 outcome, 70% complete)",
     lambda T: mdes_did(BASE.with_(completeness=0.7), T, B=6, n_comparison=150), SERIES[7], "G  DiD, 6-yr baseline"),
    ("H. As G, plus proxy covariates absorbing half the year-to-year shocks (illustrative)",
     lambda T: mdes_did(BASE.with_(completeness=0.7, k_t=0.075), T, B=6, n_comparison=150), SERIES[2], "H  DiD + proxy covariates"),
]
fig, ax = plt.subplots(figsize=(10, 6.2))
plausible_band(ax)
for label, fn, col, short in designs_sec:
    y = np.array([fn(T) for T in DURATIONS])
    ax.plot(DURATIONS, y, lw=2, color=col, label=label,
            ls="--" if label.startswith(("F", "G", "H")) else "-")
    direct_label(ax, DURATIONS[-1], y[-1], short, col)
    for T, v in zip(DURATIONS, y):
        table_rows.append((label, T, v))
style_axes(ax)
ax.set_ylim(0, 0.30)
place_labels(ax)
ax.set_title("Secondary, non-randomised analyses buy power at the price of bias risk\n"
             "Difference-in-differences against all other public facilities in DHIS2 (dashed)",
             loc="left", fontsize=11.5)
ax.legend(loc="upper right", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_secondary_nonrandomised_mdes_vs_duration.png"), dpi=170)
plt.close(fig)

# ============================================================================
# Figure 3: outcome definition (denominator / event rate) under design D
# ============================================================================
outcomes = [
    ("All under-5 deaths in catchment (65/1,000 births)", "all U5 deaths",
     BASE.with_(c0=90, baseline_years=3), SERIES[0]),
    ("Deaths at 1-59 months only (34/1,000 births; CDS-sensitive)", "1-59 month deaths",
     BASE.with_(p0=0.028 * 1.2, k_b=0.30, c0=90, baseline_years=3), SERIES[1]),
    ("Deaths within 30 days of a sick-child visit (0.6% of ~3,000 visits/yr)", "30-day post-visit deaths",
     BASE.with_(p0=0.006, m=3000, k_b=0.35, c0=90, baseline_years=3), SERIES[2]),
    ("In-facility / referral under-5 deaths only (25% of deaths captured)", "facility-recorded deaths",
     BASE.with_(completeness=0.25, k_b=0.35, c0=90, baseline_years=3), SERIES[4]),
]
fig, ax = plt.subplots(figsize=(10, 6.2))
for label, short, a, col in outcomes:
    y = np.array([mdes_parallel(a, T) for T in DURATIONS])
    ax.plot(DURATIONS, y, lw=2, color=col, label=label)
    direct_label(ax, DURATIONS[-1], y[-1], short, col)
    for T, v in zip(DURATIONS, y):
        table_rows.append(("Outcome: " + label, T, v))
style_axes(ax, ylabel="Minimum detectable relative reduction in the outcome (80% power)")
ax.set_ylim(0, 0.45)
place_labels(ax)
ax.set_title("Outcome choice matters as much as duration\n"
             "Design D (30 vs 90 facilities, 3-yr baseline adjustment) with different primary outcomes",
             loc="left", fontsize=11.5)
ax.legend(loc="upper right", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig3_outcome_definition_mdes_vs_duration.png"), dpi=170)
plt.close(fig)

# ============================================================================
# Figure 4: sensitivity to the two least-known inputs at a 3-year trial
# ============================================================================
T_REF = 3.0
kbs = np.linspace(0.10, 0.45, 15)
fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
for ax, (design_label, base_a) in zip(axes, [
        ("A. 30 vs 30, no baseline", BASE),
        ("D. 30 vs 90, 3-yr baseline adjustment", BASE.with_(c0=90, baseline_years=3))]):
    for m, col in zip([150, 250, 400], [SERIES[0], SERIES[1], SERIES[2]]):
        y = [mdes_parallel(base_a.with_(k_b=k, m=m), T_REF) for k in kbs]
        ax.plot(kbs, y, lw=2, color=col, label=f"{m} births / facility / yr")
        direct_label(ax, kbs[-1], y[-1], f"{m}", col)
    ax.axvline(BASE.k_b, color=MUTED, lw=1, ls=":")
    ax.text(BASE.k_b + 0.005, 0.02, "assumed k_b = 0.25", fontsize=8, color=MUTED)
    ax.axhspan(*PLAUSIBLE, color=INK, alpha=0.06, lw=0)
    ax.set_title(design_label, loc="left", fontsize=10.5)
    ax.set_xlabel("Between-facility CV of true under-5 mortality (k_b)")
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.set_xlim(kbs[0], kbs[-1] + 0.05)
    ax.grid(axis="x", visible=False)
axes[0].set_ylabel("Minimum detectable reduction at 3 years (80% power)")
axes[0].set_ylim(0, 0.40)
for ax in axes:
    place_labels(ax)
axes[0].legend(loc="upper left", fontsize=8.5, title="Catchment size", title_fontsize=8.5)
fig.suptitle("Sensitivity to the two inputs we know least well: facility heterogeneity and catchment size",
             x=0.01, ha="left", fontsize=11.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig4_sensitivity_kb_catchment.png"), dpi=170)
plt.close(fig)

# ============================================================================
# Figure 5: 25 vs 30 intervention facilities, power curves at 3 years
# ============================================================================
effects = np.linspace(0.0, 0.40, 81)
fig, ax = plt.subplots(figsize=(10, 5.8))
curves = [
    ("A. 30 vs 30, no baseline", lambda d: power_parallel(d, BASE, T_REF), SERIES[0], "-"),
    ("A'. 25 vs 25, no baseline", lambda d: power_parallel(d, BASE.with_(c1=25, c0=25), T_REF), SERIES[0], "--"),
    ("D. 30 vs 90, 3-yr baseline", lambda d: power_parallel(d, BASE.with_(c0=90, baseline_years=3), T_REF), SERIES[3], "-"),
    ("D'. 25 vs 90, 3-yr baseline", lambda d: power_parallel(d, BASE.with_(c1=25, c0=90, baseline_years=3), T_REF), SERIES[3], "--"),
    ("G. DiD vs 150 facilities, 6-yr baseline", lambda d: power_did(d, BASE.with_(completeness=0.7), T_REF, 6, 150), SERIES[7], "-"),
]
for label, fn, col, ls in curves:
    y = [fn(d) for d in effects]
    ax.plot(effects, y, lw=2, color=col, ls=ls, label=label)
ax.axhline(0.8, color=MUTED, lw=1, ls=":")
ax.text(0.005, 0.81, "80% power", fontsize=8, color=MUTED)
ax.axvspan(*PLAUSIBLE, color=INK, alpha=0.06, lw=0)
ax.text(PLAUSIBLE[0] + 0.002, 0.03, "plausible\ntrue effect", fontsize=8, color=INK2)
ax.set_xlabel("True relative reduction in under-5 mortality")
ax.set_ylabel("Power (two-sided alpha = 0.05)")
ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
ax.set_xlim(0, 0.40); ax.set_ylim(0, 1.0)
ax.grid(axis="x", visible=False)
ax.set_title("Power curves for a 3-year trial: 25 vs 30 intervention facilities (dashed = 25)",
             loc="left", fontsize=11.5)
ax.legend(loc="lower right", fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig5_power_curves_3yr_25_vs_30.png"), dpi=170)
plt.close(fig)

# ============================================================================
# Tables
# ============================================================================
with open(os.path.join(RES, "mdes_by_design_and_duration.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["design", "follow_up_years", "mdes_relative_reduction"])
    for r in table_rows:
        w.writerow([r[0], f"{r[1]:.2f}", f"{r[2]:.4f}"])

# power at selected effect sizes / durations for the headline designs
headline = {
    "A. 30 vs 30, no baseline": lambda d, T: power_parallel(d, BASE, T),
    "B. 30 vs 90, no baseline": lambda d, T: power_parallel(d, BASE.with_(c0=90), T),
    "C. 30 vs 30, 3-yr baseline": lambda d, T: power_parallel(d, BASE.with_(baseline_years=3), T),
    "D. 30 vs 90, 3-yr baseline": lambda d, T: power_parallel(d, BASE.with_(c0=90, baseline_years=3), T),
    "E. Stepped wedge 30, 5 steps": lambda d, T: power_stepped_wedge(d, BASE, T, 5),
    "F. DiD vs 150, 3-yr baseline": lambda d, T: power_did(d, BASE.with_(completeness=0.7), T, 3, 150),
    "G. DiD vs 150, 6-yr baseline": lambda d, T: power_did(d, BASE.with_(completeness=0.7), T, 6, 150),
}
with open(os.path.join(RES, "power_by_design_effect_duration.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["design", "follow_up_years", "true_relative_reduction", "power"])
    for name, fn in headline.items():
        for T in [1, 2, 3, 4, 5]:
            for d in [0.10, 0.15, 0.20, 0.25, 0.30]:
                w.writerow([name, T, d, f"{fn(d, T):.3f}"])

# markdown summary table for the README
lines = ["| Design | 1 yr | 2 yr | 3 yr | 4 yr | 5 yr |", "|---|---|---|---|---|---|"]
mdes_fns = {
    "A. Parallel CRT, 30 vs 30, no baseline": lambda T: mdes_parallel(BASE, T),
    "B. Parallel CRT, 30 vs 90, no baseline": lambda T: mdes_parallel(BASE.with_(c0=90), T),
    "C. Parallel CRT, 30 vs 30, 3-yr baseline adj.": lambda T: mdes_parallel(BASE.with_(baseline_years=3), T),
    "D. Parallel CRT, 30 vs 90, 3-yr baseline adj.": lambda T: mdes_parallel(BASE.with_(c0=90, baseline_years=3), T),
    "E. Stepped wedge, 30 facilities, 5 steps": lambda T: mdes_stepped_wedge(BASE, T, 5),
    "F. DiD vs 150 DHIS2 facilities, 3-yr baseline": lambda T: mdes_did(BASE.with_(completeness=0.7), T, 3, 150),
    "G. DiD vs 150 DHIS2 facilities, 6-yr baseline": lambda T: mdes_did(BASE.with_(completeness=0.7), T, 6, 150),
}
for name, fn in mdes_fns.items():
    lines.append("| " + name + " | " + " | ".join(f"{fn(T):.0%}" for T in [1, 2, 3, 4, 5]) + " |")
with open(os.path.join(RES, "mdes_summary_table.md"), "w") as f:
    f.write("\n".join(lines) + "\n")

print("\n".join(lines))
print()
print("Expected deaths, 3 yrs:", expected_events(BASE, 3))
print("Residual k_b^2 after 3-yr baseline adj:", round(residual_kb2(BASE.with_(baseline_years=3)), 4),
      "vs unadjusted", BASE.k_b ** 2)
print("Power for a true 15% reduction, design D, at 3 / 5 yrs:",
      round(power_parallel(0.15, BASE.with_(c0=90, baseline_years=3), 3), 2),
      round(power_parallel(0.15, BASE.with_(c0=90, baseline_years=3), 5), 2))
