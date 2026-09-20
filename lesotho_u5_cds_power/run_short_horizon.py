"""Short-horizon view: what is detectable between 6 and 18 months of follow-up?

Stacks the power-improving design choices cumulatively, adds an intervention
ramp-up (the tool is not at full fidelity on day one), and compares outcome
definitions against the effect size each could plausibly show.

Run:  python3 run_short_horizon.py   ->  figures/fig6_*, results/short_horizon_*.csv
"""
from __future__ import annotations

import csv
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from model import Assumptions, mdes_parallel, power_parallel
from plotting import (SERIES, INK, INK2, MUTED, AXIS, SURFACE,
                      direct_label, place_labels)

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, RES = os.path.join(HERE, "figures"), os.path.join(HERE, "results")

BASE = Assumptions()
MONTHS = np.arange(6, 18.01, 0.5)
YEARS = MONTHS / 12.0
RAMP_MONTHS = 3.0   # linear ramp to full fidelity over the first 3 months


def ramp_factor(T_years: float, ramp_years: float = RAMP_MONTHS / 12) -> float:
    """Average fraction of full effect over follow-up with a linear ramp-up."""
    if T_years <= ramp_years:
        return 0.5 * T_years / ramp_years
    return 1.0 - 0.5 * ramp_years / T_years


def mdes_with_ramp(a: Assumptions, T: float) -> float:
    return mdes_parallel(a, T) / ramp_factor(T)


def power_with_ramp(d: float, a: Assumptions, T: float) -> float:
    return power_parallel(d * ramp_factor(T), a, T)


# ---------------------------------------------------------------------------
# Cumulative stack of design choices, all-cause U5 mortality outcome
# ---------------------------------------------------------------------------
STACK = [
    ("0. Starting point: 30 vs 30, 250 births/facility, no baseline", BASE, SERIES[0]),
    ("1. + 90 control facilities", BASE.with_(c0=90), SERIES[1]),
    ("2. + adjust for 3-yr DHIS2 baseline", BASE.with_(c0=90, baseline_years=3), SERIES[2]),
    ("3. + large catchments (400 births/facility/yr)", BASE.with_(c0=90, baseline_years=3, m=400), SERIES[3]),
    ("4. + all ~150 other public facilities as controls", BASE.with_(c0=150, baseline_years=3, m=400), SERIES[4]),
]
FULL = STACK[-1][1]

# ---------------------------------------------------------------------------
# Outcome alternatives under the full stack, with the effect each could
# plausibly show (all-U5 10-15% from IMCI trials; 1-59 m scaled by the
# post-neonatal share of deaths, 28/54; post-visit mortality is a guess)
# ---------------------------------------------------------------------------
OUTCOMES = [
    ("All U5 deaths in catchment (65/1,000 births)", "all U5 deaths",
     FULL, (0.10, 0.15), SERIES[0]),
    ("Deaths at 1-59 months (34/1,000 births)", "1-59 month deaths",
     FULL.with_(p0=0.028 * 1.2, k_b=0.30), (0.19, 0.29), SERIES[1]),
    ("30-day post-visit deaths, health centres (0.6% of 3,000 visits/yr)", "post-visit deaths, HCs",
     FULL.with_(p0=0.006, m=3000, k_b=0.35), (0.20, 0.30), SERIES[2]),
    ("30-day post-visit deaths, HCs + hospital OPDs (6,000 visits/yr)", "post-visit deaths, HCs + hospitals",
     FULL.with_(p0=0.006, m=6000, k_b=0.35), (0.20, 0.30), SERIES[6]),
]

fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))

# left panel
ax = axes[0]
ax.axhspan(0.10, 0.15, color=INK, alpha=0.06, lw=0)
ax.text(6.1, 0.103, "Plausible true effect on all-cause U5 mortality: 10-15%", fontsize=8, color=INK2, va="bottom")
rows = []
for label, a, col in STACK:
    y = np.array([mdes_parallel(a, T) for T in YEARS])
    ax.plot(MONTHS, y, lw=2, color=col, label=label)
    direct_label(ax, MONTHS[-1], y[-1], label.split(". ")[0], col)
    rows += [(label, m, v) for m, v in zip(MONTHS, y)]
y_ramp = np.array([mdes_with_ramp(FULL, T) for T in YEARS])
ax.plot(MONTHS, y_ramp, lw=2, ls="--", color=SERIES[4], label="4 with a 3-month ramp-up to full fidelity")
direct_label(ax, MONTHS[-1], y_ramp[-1], "4, ramp-up", SERIES[4])
rows += [("4 with 3-month ramp-up", m, v) for m, v in zip(MONTHS, y_ramp)]
ax.set_ylim(0, 0.45)
ax.set_xlim(6, 20.5)
ax.set_xticks([6, 9, 12, 15, 18])
place_labels(ax)
ax.axvline(12, color=MUTED, lw=1, ls=":")
ax.set_xlabel("Months of intervention follow-up")
ax.set_ylabel("Minimum detectable reduction in under-5 mortality (80% power)")
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
ax.grid(axis="x", visible=False)
ax.set_title("Each design choice, stacked (all-cause under-5 mortality)", loc="left", fontsize=10.5)
ax.legend(loc="upper right", fontsize=8)

# right panel
ax = axes[1]
bar_x = [22.6, 23.2, 23.8, 24.4]
for (label, short, a, plaus, col), bx in zip(OUTCOMES, bar_x):
    y = np.array([mdes_with_ramp(a, T) for T in YEARS])
    ax.plot(MONTHS, y, lw=2, color=col, label=f"{label}; plausible effect {plaus[0]:.0%}-{plaus[1]:.0%}")
    direct_label(ax, MONTHS[-1], y[-1], short, col)
    # plausible-effect range: a short bar per outcome, in the series colour
    ax.plot([bx, bx], plaus, lw=5, color=col, alpha=0.45, solid_capstyle="butt")
    rows += [("Outcome (ramp-up): " + label, m, v) for m, v in zip(MONTHS, y)]
ax.text(23.5, 0.315, "plausible true effect,\nby outcome", fontsize=7.5, color=INK2, va="bottom", ha="center")
ax.set_ylim(0, 0.45)
ax.set_xlim(6, 25)
ax.set_xticks([6, 9, 12, 15, 18])
place_labels(ax)
ax.axvline(12, color=MUTED, lw=1, ls=":")
ax.set_xlabel("Months of intervention follow-up")
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
ax.grid(axis="x", visible=False)
ax.set_title("Outcome choice under design 4 with ramp-up", loc="left", fontsize=10.5)
ax.legend(loc="upper left", fontsize=7.5)

fig.suptitle("What can a 30-facility CDS trial in Lesotho detect within 6-18 months?",
             x=0.01, ha="left", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig6_short_horizon_6_to_18_months.png"), dpi=170)
plt.close(fig)

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
with open(os.path.join(RES, "short_horizon_mdes.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["design_or_outcome", "months", "mdes_relative_reduction"])
    for r in rows:
        w.writerow([r[0], f"{r[1]:.1f}", f"{r[2]:.4f}"])

print("MDES (80% power) by months of follow-up, with 3-month ramp-up, full design (30 vs 150, baseline adj., 400 births):")
print("| Outcome | 6 mo | 9 mo | 12 mo | 15 mo | 18 mo | plausible effect |")
print("|---|---|---|---|---|---|---|")
for label, short, a, plaus, col in OUTCOMES:
    vals = " | ".join(f"{mdes_with_ramp(a, m / 12):.0%}" for m in (6, 9, 12, 15, 18))
    print(f"| {label} | {vals} | {plaus[0]:.0%}-{plaus[1]:.0%} |")
print()
print("Power at 12 months (ramp-up included) for true effects, by outcome:")
print("| Outcome | 15% | 20% | 25% | 30% |")
print("|---|---|---|---|---|")
for label, short, a, plaus, col in OUTCOMES:
    vals = " | ".join(f"{power_with_ramp(d, a, 1.0):.0%}" for d in (0.15, 0.20, 0.25, 0.30))
    print(f"| {short} | {vals} |")
print()
print("Stack at 12 months, no ramp:", [f"{mdes_parallel(a, 1.0):.1%}" for _, a, _ in STACK],
      "| with ramp (full):", f"{mdes_with_ramp(FULL, 1.0):.1%}")
print("Full design, one-sided alpha 0.05 at 12 mo w/ ramp:", f"{mdes_with_ramp(FULL.with_(alpha=0.10), 1.0):.1%}")
print("Full design, 25 intervention facilities at 12 mo w/ ramp:", f"{mdes_with_ramp(FULL.with_(c1=25), 1.0):.1%}")
