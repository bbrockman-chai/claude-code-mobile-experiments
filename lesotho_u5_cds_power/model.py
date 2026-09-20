"""
Power / minimum-detectable-effect model for a facility-level cluster trial of
clinical decision support (CDS) aimed at under-5 mortality in Lesotho.

All designs are analysed at the cluster (facility-catchment) level using a
variance-components approximation for rare event rates:

    deaths_ij ~ Poisson( m * p * u_i * v_ij )

    m     births (or child-years / visits) per facility per year
    p     underlying event probability
    u_i   persistent facility effect,  mean 1, CV = k_b  (between-facility)
    v_ij  transient facility-year shock, mean 1, CV = k_t (year-to-year)

Variance of a facility's observed rate over T years of follow-up:

    Var(r_i) ~= p/(m T) + p^2 * ( k_b^2 + k_t^2 / T )

The three pieces are (1) Poisson sampling noise, which time buys down,
(2) persistent between-facility heterogeneity, which time does NOT buy down
and which baseline adjustment / matching can shrink, and (3) transient shocks,
which time buys down slowly.

The model is deliberately simple and transparent; every parameter is an
explicit assumption listed in README.md so it can be replaced with better data.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy import stats, optimize


# ----------------------------------------------------------------------------
# Assumptions container
# ----------------------------------------------------------------------------
@dataclass
class Assumptions:
    # outcome
    p0: float = 0.054 * 1.2          # baseline event prob per birth (U5MR x high-burden uplift)
    m: float = 250.0                 # births per facility catchment per year
    completeness: float = 1.0        # fraction of true events captured by the outcome system
    # heterogeneity
    k_b: float = 0.25                # persistent between-facility CV of the true rate
    k_t: float = 0.15                # transient facility-year CV of the true rate
    # design
    c1: int = 30                     # intervention facilities
    c0: int = 30                     # control facilities
    baseline_years: float = 0.0      # years of pre-period data used for adjustment (0 = none)
    baseline_completeness: float = 0.7  # completeness of the admin (DHIS2) baseline data
    # inference
    alpha: float = 0.05
    power: float = 0.80

    def with_(self, **kw) -> "Assumptions":
        return replace(self, **kw)


# ----------------------------------------------------------------------------
# Building blocks
# ----------------------------------------------------------------------------
def residual_kb2(a: Assumptions) -> float:
    """Persistent between-facility variance (as CV^2) left after regressing on a
    baseline rate observed over `baseline_years` of admin data.

    R^2 between the observed baseline rate and the true persistent factor u_i is
    the reliability of the baseline as a measurement of u_i.
    """
    if a.baseline_years <= 0:
        return a.k_b ** 2
    B = a.baseline_years
    p_obs = a.p0 * a.baseline_completeness
    noise = a.k_t ** 2 / B + 1.0 / (p_obs * a.m * B)   # CV^2 of noise in baseline rate
    r2 = a.k_b ** 2 / (a.k_b ** 2 + noise)
    return a.k_b ** 2 * (1.0 - r2)


def rate_variance(p: float, a: Assumptions, T: float, kb2: float) -> float:
    """Variance of one facility's observed event rate over T years of follow-up,
    on the scale of the *true* rate (i.e. divided by completeness)."""
    p_obs = p * a.completeness
    poisson = p_obs / (a.m * T) / a.completeness ** 2
    return poisson + p ** 2 * (kb2 + a.k_t ** 2 / T)


def _power_from_se(delta_abs: float, se: float, df: float, alpha: float) -> float:
    tcrit = stats.t.ppf(1 - alpha / 2, df)
    z = delta_abs / se
    return (1 - stats.norm.cdf(tcrit - z)) + stats.norm.cdf(-tcrit - z)


def _solve_mdes(power_fn, a: Assumptions) -> float:
    """Find relative reduction d in (0, 0.95) with power == target."""
    f = lambda d: power_fn(d) - a.power
    lo, hi = 1e-4, 0.95
    if f(hi) < 0:
        return np.nan  # not detectable even at 95% reduction
    return optimize.brentq(f, lo, hi, xtol=1e-5)


# ----------------------------------------------------------------------------
# Design 1: parallel cluster RCT (optionally unequal allocation, optionally
#           adjusted for a multi-year administrative baseline)
# ----------------------------------------------------------------------------
def power_parallel(d: float, a: Assumptions, T: float) -> float:
    kb2 = residual_kb2(a)
    p1 = a.p0 * (1 - d)
    var = rate_variance(a.p0, a, T, kb2) / a.c0 + rate_variance(p1, a, T, kb2) / a.c1
    df = a.c0 + a.c1 - 2 - (1 if a.baseline_years > 0 else 0)
    return _power_from_se(a.p0 - p1, np.sqrt(var), df, a.alpha)


def mdes_parallel(a: Assumptions, T: float) -> float:
    return _solve_mdes(lambda d: power_parallel(d, a, T), a)


# ----------------------------------------------------------------------------
# Design 2: stepped wedge across the intervention facilities only
#           (Hussey & Hughes 2007 variance, cluster-period means)
# ----------------------------------------------------------------------------
def sw_variance_factor(I: int, steps: int, period_years: float, a: Assumptions,
                       p: float) -> float:
    """Var(theta_hat) for a stepped wedge with I clusters, `steps` cross-over
    steps (so steps+1 periods incl. an all-control baseline period), each
    period `period_years` long."""
    Tp = steps + 1
    X = np.zeros((I, Tp))
    groups = np.array_split(np.arange(I), steps)
    for s, g in enumerate(groups):
        X[np.ix_(g, np.arange(s + 1, Tp))] = 1.0
    U = X.sum()
    W = (X.sum(axis=0) ** 2).sum()
    V = (X.sum(axis=1) ** 2).sum()
    m_period = a.m * period_years
    p_obs = p * a.completeness
    # k_t is defined per facility-year, so a period of length L years has
    # transient shock CV^2 = k_t^2 / L
    sig_e2 = p_obs / m_period / a.completeness ** 2 + a.k_t ** 2 * p ** 2 / period_years
    sig_a2 = a.k_b ** 2 * p ** 2
    num = I * sig_e2 * (sig_e2 + Tp * sig_a2)
    den = (I * U - W) * sig_e2 + (U ** 2 + I * Tp * U - Tp * W - I * V) * sig_a2
    return num / den


def power_stepped_wedge(d: float, a: Assumptions, T: float, steps: int) -> float:
    I = a.c1
    period_years = T / (steps + 1)
    p1 = a.p0 * (1 - d)
    pbar = 0.5 * (a.p0 + p1)
    var = sw_variance_factor(I, steps, period_years, a, pbar)
    df = I - steps - 1
    return _power_from_se(a.p0 - p1, np.sqrt(var), max(df, 5), a.alpha)


def mdes_stepped_wedge(a: Assumptions, T: float, steps: int = 5) -> float:
    return _solve_mdes(lambda d: power_stepped_wedge(d, a, T, steps), a)


# ----------------------------------------------------------------------------
# Design 3: non-randomised difference-in-differences against many comparison
#           facilities using administrative (DHIS2) data for both periods.
#           Persistent facility effects cancel; what remains is sampling noise
#           and transient shocks in both pre and post periods.
# ----------------------------------------------------------------------------
def power_did(d: float, a: Assumptions, T: float, B: float, n_comparison: int) -> float:
    q = a.completeness
    p1 = a.p0 * (1 - d)

    def change_var(p_post):
        pois = (a.p0 * q) / (a.m * B) / q ** 2 + (p_post * q) / (a.m * T) / q ** 2
        shock = a.k_t ** 2 * (a.p0 ** 2 / B + p_post ** 2 / T)
        return pois + shock

    var = change_var(a.p0) / n_comparison + change_var(p1) / a.c1
    df = n_comparison + a.c1 - 2
    return _power_from_se(a.p0 - p1, np.sqrt(var), df, a.alpha)


def mdes_did(a: Assumptions, T: float, B: float, n_comparison: int) -> float:
    return _solve_mdes(lambda d: power_did(d, a, T, B, n_comparison), a)


# ----------------------------------------------------------------------------
# Convenience: expected events, to sanity check
# ----------------------------------------------------------------------------
def expected_events(a: Assumptions, T: float) -> dict:
    per_fac_year = a.m * a.p0
    return {
        "events_per_facility_year": per_fac_year,
        "events_intervention_arm": per_fac_year * a.c1 * T,
        "events_control_arm": per_fac_year * a.c0 * T,
    }


# ----------------------------------------------------------------------------
# Intervention ramp-up: average fraction of full effect over follow-up T when
# the effect rises linearly from 0 to full over the first `ramp_years`.
# ----------------------------------------------------------------------------
def ramp_factor(T_years: float, ramp_years: float = 0.25) -> float:
    if T_years <= ramp_years:
        return 0.5 * T_years / ramp_years
    return 1.0 - 0.5 * ramp_years / T_years
