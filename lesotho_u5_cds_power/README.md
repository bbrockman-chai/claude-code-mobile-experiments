# Lesotho CDS under-5 mortality trial: rough power analysis

**Question.** If we can deploy a clinical decision support (CDS) tool in at most
25-30 health facilities in Lesotho, what reduction in under-5 mortality could a
trial detect, and how does that trade off against study duration?

**Short answer.** With 30 intervention facilities in high-burden catchments, a
plain 30-vs-30 cluster RCT can only detect a ~20% reduction in under-5 mortality
even after 3-5 years, which is above the 10-15% that IMCI-type interventions
have historically achieved. Three cheap design choices together roughly halve
the minimum detectable effect: (1) randomise against many more control
facilities than intervention facilities, (2) adjust for a multi-year DHIS2
mortality baseline, and (3) pick large as well as high-burden catchments. With
all three, a 3-year trial can detect ~12% (93% power for a true 15% effect). A
non-randomised difference-in-differences against all other public facilities
reaches a similar level and is a natural pre-specified secondary analysis.
Beyond ~3 years, extra duration buys very little because the binding constraint
is between-facility heterogeneity, not the number of deaths.

Everything here is driven by explicit assumptions listed below. Replace them
in `model.py` (`Assumptions` dataclass) and re-run `python3 run_analysis.py`.

---

## 1. Headline results

Minimum detectable relative reduction in under-5 mortality (80% power,
two-sided alpha 0.05, cluster-level analysis), by design and years of follow-up:

| Design | 1 yr | 2 yr | 3 yr | 4 yr | 5 yr |
|---|---|---|---|---|---|
| A. Parallel CRT, 30 vs 30, no baseline | 25% | 22% | 20% | 19% | 19% |
| B. Parallel CRT, 30 vs 90, no baseline | 20% | 17% | 16% | 15% | 15% |
| C. Parallel CRT, 30 vs 30, 3-yr baseline adj. | 22% | 18% | 16% | 15% | 14% |
| D. Parallel CRT, 30 vs 90, 3-yr baseline adj. | 18% | 14% | 12% | 12% | 11% |
| E. Stepped wedge, 30 facilities, 5 steps | 40% | 31% | 26% | 23% | 21% |
| F. DiD vs 150 DHIS2 facilities, 3-yr baseline | 20% | 16% | 15% | 14% | 13% |
| G. DiD vs 150 DHIS2 facilities, 6-yr baseline | 19% | 14% | 13% | 12% | 11% |

Power for a **true 15% reduction** at 3 years: A 52%, B 74%, C 76%, **D 93%**,
E 32%, G 92%. For a true 10% reduction none of the designs reaches 80% power
within 5 years (D reaches 71% at 5 years).

Asymptotes (infinite follow-up): design A bottoms out at ~17%, design D at ~9%.
That floor is set by residual between-facility variation, which is why the
"tricks" matter more than duration once the trial is 2-3 years long.

Going from 30 to 25 intervention facilities costs about one percentage point
of MDES under design D (12.5% -> 13.3% at 3 years).

Figures: `figures/fig1` (randomised designs), `fig2` (non-randomised secondary
analyses), `fig3` (outcome definition), `fig4` (sensitivity to heterogeneity
and catchment size), `fig5` (power curves at 3 years, 25 vs 30 facilities),
`fig6` (6-18 month horizon, section 1b), `fig7` (composite outcome, section 1c).
Tables: `results/`.

---

## 1b. Short horizon: what is doable within 12 months?

`run_short_horizon.py` zooms into 6-18 months, stacks the design choices
cumulatively and adds a 3-month linear ramp-up to full fidelity (a tool
deployed at month 0 is not changing care on day one). Figure:
`figures/fig6_short_horizon_6_to_18_months.png`.

Minimum detectable reduction in all-cause under-5 mortality at 12 months,
each choice added to the last:

| Design (cumulative) | MDES at 12 months |
|---|---|
| 0. 30 vs 30, 250 births/facility/yr, no baseline | 25% |
| 1. + 90 control facilities | 20% |
| 2. + adjust for 3-yr DHIS2 baseline | 18% |
| 3. + large catchments (400 births/facility/yr) | 15% |
| 4. + all ~150 other public facilities as randomised controls | 14% |
| 4 with a 3-month ramp-up (realistic) | 16.5% |
| 4 with ramp-up, 25 rather than 30 intervention facilities | 17.6% |
| 4 with ramp-up, one-sided test | 14.8% |

With everything stacked and a realistic ramp-up, by outcome:

| Outcome | 6 mo | 9 mo | 12 mo | 15 mo | 18 mo | plausible true effect |
|---|---|---|---|---|---|---|
| All U5 deaths in catchment | 25% | 19% | 16% | 15% | 14% | 10-15% |
| Deaths at 1-59 months | 31% | 24% | 20% | 18% | 17% | 19-29% |
| 30-day post-visit deaths, health centres (3,000 visits/yr) | 28% | 22% | 19% | 17% | 16% | 20-30% (guess) |
| 30-day post-visit deaths, HCs + hospital OPDs (6,000 visits/yr) | 23% | 18% | 15% | 14% | 13% | 20-30% (guess) |

Power at 12 months (ramp-up included) if the true effect on the chosen outcome is:

| Outcome | 15% | 20% | 25% | 30% |
|---|---|---|---|---|
| All U5 deaths | 71% | 94% | 99% | 100% |
| 1-59 month deaths | 51% | 78% | 94% | 99% |
| Post-visit deaths, HCs | 59% | 85% | 97% | 100% |
| Post-visit deaths, HCs + hospitals | 78% | 97% | 100% | 100% |

**Reading.** A 12-month result on all-cause under-5 mortality is only
convincing if the true effect is about 20% or more, which is above what IMCI
programmes have shown. Six months is not a mortality trial under any design
here (MDES 23-31%). The realistic way to get a mortality answer in 12 months
is to change the outcome to one the tool acts on directly and that has more
events per cluster: 30-day mortality after a sick-child visit, with active
follow-up (phone call or village health worker visit) in both arms, and with
hospital outpatient departments included in the clusters to double the
denominators. That outcome has a plausibly larger effect (20-30%) than
all-cause mortality and reaches ~80% power for a 15% effect and ~97% for 20%
at 12 months. The 1-59 month outcome is the right restriction conceptually
but has too few events to be quick on its own.

**Recommended 12-month package.**
1. Primary outcome: 30-day mortality after a sick-child visit (age 2-59
   months), ascertained by active follow-up independent of the CDS tool.
   Secondary: all-cause 1-59 month mortality in the catchment via village
   health worker death reporting, and process measures (correct
   classification, treatment, referral completion), which will be
   definitive within months.
2. Clusters: 30 intervention facilities chosen for burden *and* volume,
   pairing each district hospital OPD with its filter clinics; every other
   public facility in the DHIS2 panel as a randomised control (constrained
   randomisation on 3-year baseline mortality, volume and district).
3. Analysis: baseline-adjusted cluster-level comparison as primary;
   pre-specified DiD/synthetic control on the full DHIS2 panel as secondary;
   a Bayesian summary (posterior probability of >=10% reduction) reported
   alongside, which is more informative than a p-value at 12 months.
4. Timing: 3-month ramp-up, 12-month primary read-out, with a pre-specified
   option to extend to 24 months for the catchment-mortality secondary
   outcome, where the curves keep improving.
5. Seasonality: a 12-month window covers one full diarrhoea/pneumonia
   season; anything shorter must be season-matched against baseline.

Caveats specific to the short horizon: the ramp-up assumption (3 months) is a
guess and could easily be 6 for a new digital tool, which would push the
12-month numbers up by another 2-3 points; the post-visit mortality rate (0.6%)
and its between-facility CV (0.35) are guesses to be replaced by facility data;
active follow-up must achieve the same completeness in both arms.

---

## 1c. Composite outcome: neonatal deaths + 30-day post-visit deaths

`run_composite.py`, figure `figures/fig7_composite_neonatal_plus_sick_visit.png`.
The composite is "death among every child who had a facility contact": births
enrolled at delivery or first postnatal visit and followed to 28 days, plus
sick-child visits (2-59 months) followed to 30 days. One follow-up system
serves both.

Events per cluster-year (large clusters, HC + hospital OPD):

| Stream | Assumption | Events / cluster / yr |
|---|---|---|
| 30-day post-visit deaths | 6,000 visits x 0.6% | ~36 |
| Neonatal deaths | 400 births x 85% followed x 31/1,000 | ~11 |
| Composite | | ~47 |

Minimum detectable reduction at 12 months (ramp-up included, 30 vs 150,
baseline-adjusted): sick-visit only 15.2%, composite 14.2%. Without any usable
baseline both jump to 22-23%, see the caveat below.

But the composite effect is the event-weighted mean of the two component
effects, so the neonatal stream only helps if the tool moves it too.
Power at 12 months by the true effect on post-visit mortality:

| Outcome | 15% | 20% | 25% | 30% |
|---|---|---|---|---|
| Sick-visit deaths only | 79% | 97% | 100% | 100% |
| Composite, neonatal effect 0% | 61% | 87% | 98% | 100% |
| Composite, neonatal effect 10% | 78% | 95% | 99% | 100% |
| Composite, neonatal effect 20% | 90% | 98% | 100% | 100% |
| Composite, neonatal effect equal to sick-visit effect | 84% | 98% | 100% | 100% |

Break-even: the composite beats the sick-visit-only outcome when the neonatal
effect is at least ~47% of the post-visit effect. Below that it dilutes.

**Take.** Use the composite as the primary outcome only if the CDS tool has a
real newborn component (delivery and immediate newborn care, sick young infant
algorithm, postnatal danger signs) that you expect to move neonatal deaths by
at least 10%. If the tool is an IMCI-style outpatient algorithm, keep
post-visit mortality as primary and report neonatal mortality separately as a
pre-specified secondary, since including it costs up to 18 points of power
for no gain. Either way, enrolling births into the same follow-up system is
cheap and gives the neonatal answer too.

**Baseline caveat for post-visit outcomes.** DHIS2 has facility neonatal
deaths but no post-visit mortality, so the baseline adjustment for the
sick-visit stream has to work through proxies (facility U5 deaths, OPD
volume, admissions, referrals). The two curves on the left of fig 7 bracket
the answer: 15% with a baseline as good as a 3-year direct one, 23% with
none. The truth is in between and is a strong reason to run 3-6 months of
follow-up in all facilities before randomisation, which also debugs the
follow-up system.

---

## 2. Data and assumptions (all replaceable)

| Parameter | Value used | Source / reasoning |
|---|---|---|
| National U5MR | 54 per 1,000 (5 yrs before survey); 58 (10 yrs) | Lesotho DHS 2023-24, pulled via the DHS API. Trend: 113 (2004), 117 (2009), 85 (2014), 54 (2023-24). |
| Neonatal mortality | 26 per 1,000 | DHS 2023-24. So 1-59 month mortality is ~28 per 1,000, i.e. roughly half of U5 deaths are neonatal and largely outside the reach of an outpatient sick-child CDS tool. |
| District U5MR (10-yr window, 2023-24) | Thaba-Tseka 71, Leribe 67, Mokhotlong 60, Butha-Buthe 59, Berea 58, Quthing 58, Mohale's Hoek 56, Maseru 54, Mafeteng 46, Qacha's Nek 39 | DHS 2023-24 regional table. District estimates are noisy (few hundred births each). |
| High-burden uplift | x1.2 -> **65 per 1,000** used as the trial baseline | Selecting Thaba-Tseka / Leribe / Mokhotlong and the worse catchments within them. Could plausibly be 1.1-1.4. |
| Births per year, national | ~55,000-60,000 | Population ~2.3 million, crude birth rate ~25 per 1,000 (UN WPP 2024). |
| Public/CHAL facilities | ~200 health centres + ~20 hospitals | Ministry of Health / CHAL facility lists (approximate). Gives ~250 births per catchment per year on average. |
| Births per facility catchment | **250 / yr** (sensitivity 150-400) | Mountain health centres are smaller (100-200); hospitals and lowland health centres larger (400+). Catchment size is a design lever, see fig 4. |
| Expected U5 deaths per facility-year | ~16 | 250 x 0.065. Over 3 years, ~1,450 deaths per 30-facility arm. |
| Between-facility CV of true U5 mortality (k_b) | **0.25** (sensitivity 0.10-0.45) | The CV across DHS districts in 2023-24 is ~0.16 including sampling noise; facility catchments are smaller and more heterogeneous. Hayes & Moulton's rule of thumb for mortality outcomes is 0.15-0.40. This is the single most important unknown and can be estimated directly from DHIS2 facility death counts. |
| Year-to-year CV of a facility's true rate (k_t) | 0.15 | Epidemics, staffing changes, stock-outs. Guess. Raising it to 0.25 moves design D at 3 years from 12.5% to 14.3%. |
| DHIS2 completeness for U5 deaths | 70% (sensitivity 50-90%) | Community deaths are under-reported in routine data. Only matters for the baseline adjustment and the DiD designs. At 50% completeness, design G at 3 years moves from 13% to 14.5%. |
| Plausible true effect | 10-15% reduction in all-cause U5 mortality | Tanzania IMCI effectiveness study (Armstrong Schellenberg et al., Lancet 2004): 13% lower U5 mortality. Bangladesh Matlab IMCI cluster RCT (Arifeen et al., Lancet 2009): 13% reduction, not significant. Electronic CDS trials (ALMANACH, ePOCT+, e-IMCI) have measured quality of care and antibiotic use, not mortality. |
| Test | Two-sided alpha 0.05, 80% power, t-distribution with cluster-level degrees of freedom | Standard. One-sided or alpha 0.10 would shave 2-3 points off every MDES and are noted, not used. |

### How the model works

Deaths in facility i, year j ~ Poisson(m x p x u_i x v_ij), with u_i a
persistent facility factor (CV k_b) and v_ij a transient facility-year shock
(CV k_t). The variance of a facility's observed rate over T years is

    p/(mT)  +  p^2 * ( k_b^2 + k_t^2 / T )

Time buys down the first term (Poisson noise) and, slowly, the third. It does
nothing to the second, which is why every curve flattens. Baseline adjustment
shrinks k_b^2 to k_b^2 (1 - R^2), where R^2 is the reliability of the observed
baseline rate as a measure of u_i; longer and more complete baselines raise
R^2. The stepped wedge uses the Hussey & Hughes (2007) variance with the same
components. The DiD design cancels u_i entirely but pays Poisson and shock
variance twice (pre and post) and inherits DHIS2 completeness.

---

## 3. Design alternatives and the "tricks", ranked by value

1. **Unequal allocation: 30 intervention vs 60-90 control facilities.**
   Control facilities cost almost nothing if the outcome comes from routine
   data or a light surveillance system. Moving from 30 to 90 controls takes
   the 3-year MDES from 20% to 16% on its own. Diminishing returns past ~3:1.

2. **Multi-year DHIS2 baseline as a covariate (ANCOVA), and restricted or
   pair-matched randomisation on it.** Three years of routine facility death
   counts plus size and district cut residual between-facility variance by
   about 60% under the assumptions here (k_b^2 from 0.063 to 0.023). Six years
   of baseline instead of three helps only a little more (12.5% -> 11.4%).
   Covariate-constrained randomisation on baseline mortality, catchment size
   and district makes the arms balanced by construction, which also protects
   credibility with 30 clusters.

3. **Pick large catchments as well as high-burden ones.** Under design D at
   3 years, 150 vs 400 births per catchment is the difference between a 15%
   and an 11% MDES (fig 4, right). If the tool is aimed at primary care,
   pairing each hospital's outpatient department with its filter clinics as
   one cluster is a way to get big clusters without diluting the intervention.

4. **Outcome definition.** All-cause U5 mortality in the catchment is the
   cleanest but half of it is neonatal and mostly outside the tool's reach.
   Options (fig 3):
   - *1-59 month mortality* as primary: fewer events, MDES 16% at 3 years, but
     the plausible true effect on this outcome is larger (a 13% all-U5 effect
     concentrated in the post-neonatal half is a ~25% post-neonatal effect).
     This is the better-powered choice relative to the plausible effect.
   - *30-day mortality after a sick-child visit*, with active phone or village
     health worker follow-up: many more denominators, similar MDES to all-U5,
     but a much more direct test of what CDS does and cheaper to measure
     completely in both arms.
   - *Facility-recorded deaths only*: 25% capture makes this weak (MDES 19%
     at 3 years) and, worse, a CDS tool that improves record-keeping will
     increase recorded deaths in the intervention arm. Avoid as a primary
     outcome.
   - A composite (death or hospital admission for severe illness) roughly
     triples events and is worth pricing.

5. **Non-randomised difference-in-differences against the other ~150-180 public
   facilities (secondary analysis, design F/G).** Reaches about the same MDES
   as design D. Pre-specify it with: a 6-year DHIS2 baseline, facility fixed
   effects, district-by-year effects, proxy covariates that absorb
   facility-year shocks (OPD volume, immunisation coverage, ANC volume,
   stock-out days, HIV testing volumes, rainfall or malaria/diarrhoea season
   indices), placebo-in-time and pre-trend tests, and a synthetic-control
   variant as robustness. Halving the shock variance with proxies is
   illustrative (design H) and worth perhaps one point of MDES. The bias risks
   are selection of intervention sites, differential changes in reporting
   completeness caused by the tool itself, and spillover to neighbours. It
   cannot replace the randomised comparison, but it can borrow strength for
   the primary analysis if the trial's own control facilities also sit in the
   DHIS2 panel.

6. **Stepped wedge across the 30 facilities alone: not recommended for power.**
   With a rare outcome and a low ICC, the within-facility comparison is weak;
   MDES is 26% at 3 years. It only makes sense if there is no way to hold
   control facilities, and then a long pre-period from DHIS2 is essential.

7. **Other levers, smaller:** individual-level covariates (age, season) in a
   mixed model; a pre-specified interim futility look at 18-24 months; a
   Bayesian analysis with a sceptical prior from the IMCI trials so the
   result is reported as a posterior probability of benefit rather than a
   binary test; one-sided testing if the tool is judged unlikely to cause
   harm (each worth 1-3 points of MDES).

---

## 4. Things this rough model does not do

- No contamination or care-seeking across catchments (mountain catchments are
  fairly self-contained; lowland ones are not).
- No migration or in-facility referral deaths being counted at the receiving
  hospital rather than the origin catchment.
- Cluster-level analysis with a normal approximation; a mixed-model or GEE
  analysis on individual data would be a few percent more efficient.
- Assumes the outcome is measured with the same completeness in both arms.
  This is the main threat to validity for any outcome the CDS tool itself
  records. A community death surveillance system run by village health
  workers, independent of the tool, is the safest primary measurement.
- No cost. Roughly, the marginal cost of control facilities is the outcome
  measurement, so the 30-vs-90 design is only "free" if outcomes come from
  routine data or a low-cost surveillance layer.

---

## 5. Next data to fetch (in order of how much they move the answer)

1. Facility-level under-5 death counts and catchment populations from DHIS2
   for the last 6 years: gives k_b, k_t and completeness directly.
2. Facility list with catchment births or under-1 population, by district.
3. Sick-child outpatient visits per facility per year and any post-visit
   mortality follow-up data (for the visit-based outcome).
4. Village health worker death reporting coverage by district.
5. Any Lesotho-specific evidence on IMCI adherence and the share of U5 deaths
   that involved a facility contact in the preceding 30 days.

---

## Files

- `model.py`: variance-components power model (parallel CRT with unequal
  allocation and baseline adjustment, stepped wedge, DiD).
- `run_analysis.py`: produces figures 1-5 and `results/`.
- `run_short_horizon.py`: figure 6 and `results/short_horizon_mdes.csv` (6-18 months, ramp-up).
- `run_composite.py`: figure 7 and `results/composite_mdes.csv` (neonatal + post-visit composite).
- `plotting.py`: shared chart styling.
- `results/mdes_by_design_and_duration.csv`, `results/power_by_design_effect_duration.csv`,
  `results/mdes_summary_table.md`.

Run: `pip install numpy scipy matplotlib` then `python3 run_analysis.py`.
