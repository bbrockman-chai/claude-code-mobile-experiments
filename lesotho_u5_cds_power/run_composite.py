"""Composite outcome: neonatal deaths (births followed to 28 days) plus deaths
within 30 days of a sick-child visit, i.e. mortality among every child who had
a facility contact. Compares it with the sick-visit-only outcome under the full
short-horizon design (30 vs 150 facilities, large clusters, ramp-up).

Modelled as one event stream per cluster-year with expected count
lambda = E_neonatal + E_sick; the existing variance-components model applies
with m = 1 and p0 = lambda. The composite effect is the event-weighted mean of
the component effects, so adding a stream the tool barely moves dilutes power.

Run:  python3 run_composite.py  ->  figures/fig7_*, results/composite_*.csv
"""
from __future__ import annotations

import csv
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from model import Assumptions, mdes_parallel, power_parallel, ramp_factor
from plotting import SERIES, INK, INK2, MUTED, direct_label, place_labels

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, RES = os.path.join(HERE, "figures"), os.path.join(HERE, "results")

# --- event streams per cluster-year (large clusters: HC + hospital OPD) ------
VISITS = 6000            # sick-child visits (2-59 months) per cluster per year
P_VISIT_DEATH = 0.006    # 30-day mortality after a sick-child visit (guess)
BIRTHS = 400             # births per cluster per year
FOLLOWED = 0.85          # share of births enrolled and followed to 28 days
NMR = 0.026 * 1.2        # DHS 2023-24 neonatal mortality x high-burden uplift
E_SICK = VISITS * P_VISIT_DEATH          # ~36 deaths / cluster-yr
E_NEO = BIRTHS * FOLLOWED * NMR          # ~10.6 deaths / cluster-yr
E_COMP = E_SICK + E_NEO

MONTHS = np.arange(6, 18.01, 0.5)
YEARS = MONTHS / 12


def design(lam: float, k_b: float, baseline_years: float = 3) -> Assumptions:
    return Assumptions(p0=lam, m=1.0, k_b=k_b, k_t=0.15, c1=30, c0=150,
                       baseline_years=baseline_years, baseline_completeness=0.8)


SICK = design(E_SICK, 0.35)
COMP = design(E_COMP, 0.33)          # two partly independent streams -> slightly lower CV
SICK_NB = design(E_SICK, 0.35, 0)    # no usable baseline (DHIS2 has no post-visit deaths)
COMP_NB = design(E_COMP, 0.33, 0)


def mdes_r(a, T):
    return mdes_parallel(a, T) / ramp_factor(T)


def composite_effect(d_sick, d_neo):
    return (E_SICK * d_sick + E_NEO * d_neo) / E_COMP


# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 6))

# Left: MDES on each outcome's own scale, 6-18 months
ax = axes[0]
rows = []
for label, short, a, col, ls in [
    ("Sick-visit deaths only (~36 events/cluster/yr), 3-yr proxy baseline", "sick only", SICK, SERIES[0], "-"),
    ("Composite: + neonatal deaths (~47 events/cluster/yr), 3-yr baseline", "composite", COMP, SERIES[1], "-"),
    ("Sick-visit only, no usable baseline", "sick only, no baseline", SICK_NB, SERIES[0], "--"),
    ("Composite, no usable baseline", "composite, no baseline", COMP_NB, SERIES[1], "--"),
]:
    y = np.array([mdes_r(a, T) for T in YEARS])
    ax.plot(MONTHS, y, lw=2, color=col, ls=ls, label=label)
    direct_label(ax, MONTHS[-1], y[-1], short, col)
    rows += [(label, m, v) for m, v in zip(MONTHS, y)]
ax.set_ylim(0, 0.40)
ax.set_xlim(6, 21.5)
ax.set_xticks([6, 9, 12, 15, 18])
place_labels(ax)
ax.axvline(12, color=MUTED, lw=1, ls=":")
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
ax.grid(axis="x", visible=False)
ax.set_xlabel("Months of intervention follow-up")
ax.set_ylabel("Minimum detectable reduction in the outcome (80% power, with ramp-up)")
ax.set_title("Detectable effect on each outcome's own scale", loc="left", fontsize=10.5)
ax.legend(loc="upper right", fontsize=8)

# Right: power at 12 months vs the true effect on sick-visit deaths
ax = axes[1]
T = 1.0
ds = np.linspace(0.0, 0.35, 71)
ax.plot(ds, [power_parallel(d * ramp_factor(T), SICK, T) for d in ds], lw=2.6, color=SERIES[0],
        label="Sick-visit deaths only")
neo_scen = [(0.0, SERIES[7], "composite, neonatal effect 0%"),
            (0.10, SERIES[3], "composite, neonatal effect 10%"),
            (0.20, SERIES[2], "composite, neonatal effect 20%"),
            (None, SERIES[1], "composite, neonatal effect = sick-visit effect")]
prow = []
for dn, col, label in neo_scen:
    y = [power_parallel(composite_effect(d, d if dn is None else dn) * ramp_factor(T), COMP, T) for d in ds]
    ax.plot(ds, y, lw=2, color=col, ls="--", label=label.capitalize())
    prow.append((label, y))
ax.axhline(0.8, color=MUTED, lw=1, ls=":")
ax.text(0.30, 0.76, "80% power", fontsize=8, color=MUTED)
ax.axvspan(0.20, 0.30, color=INK, alpha=0.06, lw=0)
ax.text(0.202, 0.03, "plausible effect on\npost-visit mortality", fontsize=8, color=INK2)
ax.set_xlim(0, 0.35); ax.set_ylim(0, 1)
ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
ax.grid(axis="x", visible=False)
ax.set_xlabel("True reduction in 30-day post-visit mortality")
ax.set_ylabel("Power at 12 months (two-sided alpha 0.05, ramp-up included)")
ax.set_title("Composite helps only if the tool also moves neonatal deaths", loc="left", fontsize=10.5)
ax.legend(loc="upper left", fontsize=8, bbox_to_anchor=(0.0, 0.97))

fig.suptitle("Composite outcome: neonatal deaths + deaths within 30 days of a sick-child visit "
             "(30 vs 150 facilities, large clusters)", x=0.01, ha="left", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig7_composite_neonatal_plus_sick_visit.png"), dpi=170)
plt.close(fig)

with open(os.path.join(RES, "composite_mdes.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["outcome", "months", "mdes_relative_reduction"])
    for r in rows:
        w.writerow([r[0], f"{r[1]:.1f}", f"{r[2]:.4f}"])

# ---- console summary ------------------------------------------------------
be = (np.sqrt(E_SICK * E_COMP) - E_SICK) / E_NEO
print(f"Events per cluster-year: sick-visit {E_SICK:.1f}, neonatal {E_NEO:.1f}, composite {E_COMP:.1f}")
print(f"MDES at 12 mo (ramp-up): sick-only {mdes_r(SICK, 1):.1%}, composite {mdes_r(COMP, 1):.1%}; "
      f"without baseline: {mdes_r(SICK_NB, 1):.1%} / {mdes_r(COMP_NB, 1):.1%}")
print(f"Break-even: composite beats sick-only when neonatal effect >= {be:.0%} of the sick-visit effect")
print("Power at 12 mo | sick effect: 15% | 20% | 25% | 30%")
print("sick-only      |", " | ".join(f"{power_parallel(d * ramp_factor(1), SICK, 1):.0%}" for d in (0.15, 0.2, 0.25, 0.3)))
for dn, _, label in neo_scen:
    print(f"{label:45s}|", " | ".join(
        f"{power_parallel(composite_effect(d, d if dn is None else dn) * ramp_factor(1), COMP, 1):.0%}"
        for d in (0.15, 0.2, 0.25, 0.3)))
