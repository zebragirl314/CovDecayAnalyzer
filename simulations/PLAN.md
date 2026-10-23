# simulations/PLAN.md: Track A2 simulation study

Status: revised at the start of A2-1 (2026-10-22) to reflect Carly's
resolved decisions on Q1 to Q10 (section 9, "Resolved"). The scope is a
parameter space spanning plausible dsDNA-virus scenarios with localized
replication origins; it is not HCMV-calibrated. Section 10 M1 is a
post-A2-1 decision (see `NOTES.md`).

Scope: the simulation study described in the "Simulation study" section of
`NOTES.md`. The package work (`PLAN.md`), the methods paper (Track B), and
pipeline integration (Track C) are out of scope here.

Conventions used below: "decay rate" `s` is positive when coverage falls with
distance from the anchor. The LMM slope coefficient is `-s`. Effect sizes are
reported as a difference in decay rate `Delta` (per bp), so a positive
`Delta` means the treated group decays faster. Statements marked **[to
verify]** are predictions from reading the code, not measured results; the
pilot (section 3, Tier 0) is where they get checked.

---

## 1. Scope: the claims the study must support

Each claim is stated so it can fail. The estimand and the acceptance rule are
both fixed before any Tier 2 or Tier 3 run (see section 9, Q7).

**C1. Type I error of the primary test (omnibus).** Under the primary
generative model with equal decay rates across groups, the Kenward-Roger F
test on `distance_from_anchor:group` rejects at rate close to nominal
`alpha`, for `alpha` in {0.01, 0.05, 0.10}, across samples per group
`n = 2..10` and groups `G = 2, 3, 4`. Estimand: `P(p <= alpha | H0)`. We
store p-values, not rejection flags, so any `alpha` and a uniformity check
(ECDF, QQ plot) come for free.

**C1b. Type I error of the post-hoc contrasts.** For `G >= 3`, under the
global null and under *partial* nulls (some groups share a rate, one differs),
the family-wise rate of any false Tukey-adjusted pairwise rejection among the
truly equal pairs, both ungated and gated on the omnibus p < 0.05 (which is
how `test_decay_slopes` documents using them).

**C2. Power.** `P(p <= 0.05)` as a function of `n`, `G`, effect size, and
nuisance variance, over a range of effect sizes spanning "undetectable" to
"obvious", including the smallest `n` that reaches 80% power per scenario.

**C3. Head-to-head against three benchmark arms.** On identical simulated
datasets (paired by construction), compare Type I error and power of the
primary LMM against three benchmarks, omnibus and pairwise: (i) KW (chi-square
p-value) + Dunn-Holm, (ii) KW (exact permutation p-value) + Dunn-Holm, and
(iii) a two-stage parametric arm (per-sample slope, then Welch's t for `G = 2`
or one-way ANOVA for `G >= 3`). Dunn-Holm is computed both gated on KW `p <
0.05` and ungated. Paired differences get Monte Carlo CIs (McNemar style), not
two independent rates. A random-slope-only LMM comparator
(`(0 + distance_scaled | sample)`) is a conditional fourth arm, triggered only
if Tier 0 confirms M1. Benchmarks are implemented in A2-2.

**C4. Sensitivity.** How much C1 and C2 change under (a) residual
autocorrelation and mean-curve misspecification, (b) zero-coverage dropouts
(parked item 8), (c) heteroscedastic random slopes or residuals across
groups. Reported as change relative to the baseline cell, not as new
standalone claims.

**C5. Power surface.** A fitted surface for the sample size calculator, with a
stated domain of validity and a held-out accuracy check (target: held-out
absolute error on power below 0.03 inside the domain).

**C6. Anchor-specificity test (secondary). In scope at reduced size, Tier 4
(Q3). Not implemented in A2-1.** Type I error and power of the `(k+1)/(n+1)` empirical p-value from
`test_decay_significance()`, under the specificity null ("a real group effect
exists, but it is not localized at the observed anchor") and under a localized
effect. This is the expensive claim: each replicate costs `n_random + 1` LMM
fits.

What the study does *not* claim: validity on any specific organism's real data
(HCMV included; that is Track B's empirical section), validity for circular
genomes, or anything about reads-per-million quantification. The study
characterizes behavior across a parameter space; HCMV is one application point
within it, not its center.

---

## 2. Generative model

### 2.1 Design principle: simulate before normalization, reuse the package

The single most consequential choice. The package does not analyze raw
coverage; it analyzes `log(mean_cov / median(mean_cov) + 1e-5)` per sample
(`normalize_coverage`). Median normalization removes sample depth and, as a
side effect, pins each sample's intercept to its slope (see M1). A generator
that draws log-scale values directly and skips normalization would give the
wrong random-effects structure and therefore the wrong Type I and power.

So the generator:

1. Draws per-sample **mean coverage on a 1 kb "atom" grid** (before any
   binning or normalization).
2. Aggregates atoms to the analysis bin size with the same arithmetic as
   `bin_coverage` (mean over atoms; `bin_mid = (start + end) / 2`; final bin
   truncated at `genome_length`).
3. Calls the package's own `normalize_coverage` (reached with `:::`) to get
   `norm_cov` and `log_norm_cov`.
4. Attaches `sample` and `group`, and hands the table to
   `test_decay_slopes()` / `test_decay_significance()` unchanged.

Skipping per-base files saves roughly a thousand-fold in memory and time;
steps 2 and 3 are the only parts of `load_all_samples` that affect the
statistics. A one-off check (in `scripts/check_generator.R`) confirms that
writing a simulated sample as per-base `Exported.coverage.tab`, then loading
it with `load_all_samples`, reproduces the atom-path `log_norm_cov` to
numerical tolerance.

Generating on a common 1 kb atom grid makes bin size a *paired* factor: the
same simulated datasets can be analyzed at 1 kb and 5 kb at the cost of one
extra fit, not a new simulation.

### 2.2 Baseline model

For sample `i` in group `g(i)`, atom `a` at midpoint `x_a`, with true anchor
`a*`:

```
mean_cov_ia = D_i * exp( c_i - s_i * |x_a - a*| + e_ia )

s_i  = s_ref * (1 + Delta_std * I[g(i) = treated] + u_i),  u_i ~ N(0, cv^2)   (between-sample slope variation)
D_i  ~ LogNormal(log 200, 0.4^2)                              (depth; removed by normalization)
c_i  = 0                                                      (see below)
e_ia ~ N(0, sigma_atom^2), independent across atoms and samples (within-sample noise)
```

- **Decay.** Log-linear in absolute distance, symmetric arms, matching what
  the LMM assumes and what `inst/extdata` uses (`500 * exp(-decay * |pos -
  anchor|)`). The baseline is deliberately *correctly specified*, so C1 and
  C2 are statements about the test, not about robustness. Misspecification
  lives in the sensitivity analyses.
- **Between-sample variability.** An additive random slope `u_i`. There is
  no free per-sample intercept: normalization eliminates it (M1), and the
  depth factor `D_i` is included only so the simulation exercises the
  normalization step rather than assuming it away. `tau` is the quantity
  that governs power at small `n`, because the group slope difference is
  tested against between-*sample* slope variation, not against bin noise
  (see 2.5).
- **Within-sample noise.** Independent Gaussian noise on the log scale at the
  atom level. Poisson counting noise on 1 kb means is negligible at plausible
  depths (relative SD about 0.2% at 200x), so it is off in the baseline and
  switched on only in the low-depth sensitivity (`D` median 10). One
  validation figure (appendix) compares log-normal atoms against a per-base
  negative-binomial generator at three depths to justify the approximation.
- **`s_ref` is a scale parameter (adopted after A2-1).** Slopes are
  `s_i = s_ref * (1 + Delta_std * treated + u_i)` with `u_i ~ N(0, cv^2)`, so
  `tau = s_ref * cv` and `Delta = s_ref * Delta_std` are derived, and
  `Delta_std` (relative effect) and `cv` (slope coefficient of variation) are
  the generative inputs. Changing `s_ref` rescales every slope-scale quantity
  together, instead of leaving `tau` and `Delta` fixed in per-bp units. In the
  pre-reparameterization generator `tau_eff` did not depend on `s_ref`; now it
  does (2.5). Reference values: `s_ref = 6e-6`, `cv = 1/3` (`tau = 2e-6`).
  `Delta_std` below replaces `Delta / s_ref`; the older `Delta / tau_eff` in
  2.5 and 4 is a different quantity and keeps its name `delta_std` in result
  files as `Delta / tau_eff`.
- **Group effect (one group differs).** `Delta_g` shifts the group's mean
  decay rate; group 1 is the reference. Effect size `Delta` is the difference
  between the one differing group and the shared reference rate. Exact
  parameterization (Q9): `G = 2`: groups have slopes scaled by 1 and `1 + Delta_std`.
  `G = 3`: two groups have scale 1, one has `1 + Delta_std`. `G = 4`: three
  groups have scale 1, one has `1 + Delta_std`. The differing group is the
  last one (`g = G`). Equal spacing is not used. Because the omnibus F has `G
  - 1` numerator df, a single differing group is diluted as `G` grows; power
  falling with `G` at fixed `Delta` is expected, not a defect. The same
  pattern with `Delta = 0` is the global null; partial nulls for C1b are
  patterns in which the differing group is absent from the pair under test.
- **Reference rate.** `s_ref` is a primary axis (section 4), with levels 3e-6,
  6e-6 and 1.2e-5 per bp, spanning a twofold coverage drop over roughly 230 kb
  to 58 kb. The middle level, 6e-6 (twofold over about 115 kb), is the center
  of the plausible range, not an estimate for any organism. The shipped
  example data (5e-4 vs 3e-3 per bp on a 2 kb genome) is outside the range and
  is not used for calibration.
- **Slope sign.** `s_i = s_ref + Delta_g + u_i` is not truncated, so at the
  low `s_ref` and high `tau` corners a minority of samples have negative `s_i`
  (coverage rising with distance). This keeps the random-slope distribution
  Gaussian, as the LMM assumes. The realized fraction is a function of
  `s_ref / tau` and is stored per scenario.

### 2.3 How the sensitivity analyses modify the baseline

Each sensitivity is a one-factor-at-a-time modification of the baseline,
crossed only with a small set of reference cells (section 4.3). All
sensitivity datasets are generated from the same seeds as the baseline
datasets wherever possible, so deltas are paired.

**S1. Residual autocorrelation and mean-curve misspecification** (NOTES item
8, first half). Two distinct mechanisms, because they differ in what the LMM
should do about them:

- *S1a, stationary autocorrelation.* `e_ia` becomes AR(1) along the genome
  with coefficient `phi` in {0, 0.3, 0.6, 0.9} per atom, marginal SD held at
  `sigma_atom`. Isolates the "residuals not independent" assumption.
- *S1b, systematic misfit.* The true mean curve is not log-linear: coverage
  `D_i * (f + (1 - f) * exp(-s_i * d))` with floor `f` in {0.0 (baseline),
  0.2, 0.4}, a replicating-over-unreplicated-templates shape. Residuals from
  the log-linear fit are then autocorrelated *because* of the misfit, which
  is the mechanism NOTES item 8 actually describes. A sample-specific
  variant adds random `f_i`. Under S1b the interpretation of `Delta` changes
  (it is a rate difference in the true curve, and the LMM slope is a biased
  summary of it), so S1b reports power for detecting a difference in true
  `s` and states that the LMM slope is not on that scale.

**S2. Zero-coverage dropouts** (NOTES item 8, second half). Dropouts are
generated as *runs* on the atom grid (deletions and unmappable stretches are
contiguous, not independent), by a per-sample Poisson process of run starts
with geometric run lengths, setting atom coverage to exactly zero before
binning. Because binning averages atoms, a 5 kb bin goes to exactly zero (and
`log_norm_cov` to about -11.5) only if all five atoms are zeroed; partial
overlap just lowers the bin. That makes bin size matter here, so S2 is run at
both bin sizes. Axes:

| Factor | Levels |
|---|---|
| Fraction of atoms zeroed per affected sample | 0, 0.5%, 2%, 5% |
| Mean run length | 1 kb, 5 kb |
| Placement | uniform; enriched at the distal third; group-differential (rate ratio 1 vs 3) |
| Samples affected | all; one sample only |
| Analysis variant | as-is (package default); zero bins dropped in the driver before the test |

The "group-differential" placement is the dangerous one for Type I error
(dropouts correlated with group mimic an effect); the "distal" placement is
the dangerous one for slope bias (highest leverage). The drop-zero-bins
variant is implemented in the simulation driver only, as a diagnostic for what
a recommendation could look like. It does not change the package (see M5).

**S3. Heteroscedasticity across groups.** (a) Random-slope SD proportional to
the group's slope, `tau_g = cv * s_g`, instead of common `tau`; (b) residual
SD ratio between groups in {1, 2}. The LMM assumes common variance
components; with `n = 2` per group this is a plausible failure point.

**S4. Low depth.** Median `D` = 10 with the Poisson counting layer on.

**S5. Anchor misplacement.** Analysis anchor offset from the true anchor by
{0, 5, 20, 50} kb. Cheap, and it is what a user choosing OriLyt from a
reference coordinate actually risks.

### 2.4 Replicate-level quantities that are always stored

Per replicate and per analysis arm: `p_value`, `F_stat`, `df1`, `df2`,
`singular`, `conv_code`, status (`ok` or the error message), `sigma(fit)`, the
estimated random-slope SD, the realized group slope estimates, wall time, and
the scenario ID and replicate ID. The fitted variance components are stored so
users' pilot estimates can be mapped onto the surface (2.5).

### 2.5 What the generator implies for the grid (a hypothesis to test early)

The F statistic compares group slope differences against the variance of a
sample-level slope, approximately `tau^2 + sigma_bin^2 / Sxx` with `tau = s_ref * cv`, where
`Sxx = sum_j (d_j - dbar)^2` over a sample's bins. If that is right, then
`sigma`, bin size, genome length and anchor position matter for power
*only through* this one effective between-sample slope variance, and power
collapses onto `(n, G, Delta / tau_eff)`. If Tier 1 confirms the collapse,
the sample size calculator needs one user-supplied number (an estimate of
`tau_eff` from a pilot `VarCorr`) instead of five, and the grid in 4 can be
much smaller than a full factorial. If it does not collapse, the calculator is
harder and the grid stays wide. Testing this is the first job of Tier 1.

---

## 3. Replicates and tiers

Monte Carlo SE of an estimated proportion is `sqrt(p(1-p)/R)`.

| Quantity | `p` | R = 200 | R = 500 | R = 1000 | R = 2000 | R = 5000 |
|---|---|---|---|---|---|---|
| Type I at alpha = 0.05 | 0.05 | 0.015 | 0.0097 | 0.0069 | 0.0049 | 0.0031 |
| Power near 80% | 0.8 | 0.028 | 0.018 | 0.013 | 0.009 | 0.006 |
| Power near 50% | 0.5 | 0.035 | 0.022 | 0.016 | 0.011 | 0.007 |

At `R = 200` a Type I estimate of 0.07 and one of 0.03 are statistically
indistinguishable from 0.05. Telling 0.05 from 0.065 needs a CI half-width
near 0.006, i.e. `R = 5000`. Power curves are smoother in effect size than
Type I is in anything, and the surface fit pools information across cells, so
moderate `R` over many cells beats large `R` over few.

| Tier | Purpose | R per cell | Output status |
|---|---|---|---|
| 0. Pilot | Timing, failure rates, singular rates, generator checks, check the collapse hypothesis (2.5) crudely | 20 | CSVs and a short report in `results/explore/pilot/`; diagnostic, not a claim |
| 1. Exploratory | Locate interesting regions, decide which axes matter, freeze the grid | 200 | Committed as `results/explore/` |
| 2. Final Type I | C1, C1b, C3 (Type I half), C4 Type I | 5000 baseline; 2000 sensitivity | Committed |
| 3. Final power | C2, C3 (power half), C4 power | 1000 on the figure cells; 300 on the surface design | Committed |
| 4. Secondary (C6, if approved) | Anchor-specificity | 300, `n_random = 199` | Committed |

Process rule: the generative model and the analysis arms are frozen (a single
commit, referenced in every result file) before Tier 2 starts. A change to
either after that point invalidates Tier 2 and 3 results and requires
rerunning, which the `targets` dependency tracking in section 5 enforces.

### Compute estimates (every unit-time figure is a placeholder)

These are the pre-A2-1 placeholders; the A2-1 session log replaces them with
measured values. I had not timed a single fit. `R/null_distribution.R` documents one to a few
seconds per fit for about 30 samples; for the 4 to 40 samples used here I
assume roughly 0.3 s per `test_decay_slopes()` call including the `emmeans`
steps. Tier 0 replaces this with a measurement.

| Tier | Cells | Fits | CPU-hours at 0.3 s/fit | Wall time, 10 local workers |
|---|---|---|---|---|
| 1 | about 150 (see 4.2) | about 30,000 | 2.5 | 15 min |
| 2 | about 21 baseline + about 60 sensitivity, x2 bin sizes | about 400,000 | 33 | 3.5 h |
| 3 | about 105 figure cells + 400 surface points | about 225,000 | 19 | 2 h |
| 4 | about 6 cells | `6 x 300 x 200` = 360,000 | 30 | 3 h |

The local Mac (12 logical cores) can do Tiers 0, 1 and probably 2 and 3
overnight. Tier 4 and any widening of the grid are where Discovery would start
to matter; that question stays parked as instructed. The plan below keeps the
pipeline executor-agnostic so moving to HPC changes the launcher, not the
targets.

---

## 4. Parameter grid

The study maps behavior over a parameter space; it is not centered on any one
organism. Parameter values below are "centers of plausible ranges" for dsDNA
viruses with localized replication origins.

### 4.1 Axes and roles

**Design axes** (varied in the main grid):

| Axis | Values | Notes |
|---|---|---|
| Samples per group `n` | 2, 3, 4, 5, 6, 8, 10 (Tier 1 subset: 2, 3, 5, 10) | Balanced groups in the main grid. Unbalanced designs (2 vs 5; 3 vs 8, `G = 2`) are a sensitivity family. |
| Number of groups `G` | 2, 3, 4 | One-differs pattern (2.2). |
| Relative effect `Delta_std` | 0, 1/6, 1/3, 1/2, 3/4, 1 | Treated slope is `s_ref * (1 + Delta_std)`. At `s_ref = 6e-6` these are `Delta` = 0, 1, 2, 3, 4.5, 6e-6 per bp, the earlier absolute grid. |

**Geometry axes** (promoted to primary; treated as nuisance variables for
sensitivity mapping, not crossed with the design axes):

| Axis | Values | Reference value |
|---|---|---|
| Genome length `L` | 50,000; 100,000; 238,080; 500,000 bp | 238,080 |
| Reference decay rate `s_ref` (per bp) | 3e-6, 6e-6, 1.2e-5 | 6e-6 |
| Anchor position (fraction of `L`) | 0.1, 0.25, 0.5, 0.75, 0.9 | 0.5 |

`s_ref` levels give a twofold coverage drop over about 230, 115 and 58 kb. The
anchor reference of 0.5 is the center of the range (interior, two equal
arms); HCMV's OriLyt (about 0.39 L) is one application point inside the range,
not the reference. The extremes, anchor at 0 or at `L`, have one arm and are a
separate sensitivity (S6 below), not part of the five levels.

**Nuisance axes** (one-at-a-time around the reference, then crossed at reduced
size):

| Axis | Values (reference in bold) | Notes |
|---|---|---|
| Slope CV `cv` (between-sample SD = `s_ref * cv`) | 1/6, **1/3**, 2/3 | Equals `tau` = 1e-6, 2e-6, 4e-6 at `s_ref = 6e-6`. Governs power through `Delta / tau_eff` (2.5). |
| Atom-level residual SD `sigma_atom` | 0.05, **0.15**, 0.30, 0.50 | At 1 kb atoms; SD at 5 kb bins is `sigma_atom / sqrt(5)` for independent noise. Realized `sigma(fit)` is stored. |

Interpretation of `Delta = s_ref * Delta_std` at the reference (`s_ref = 6e-6`, `cv = 1/3`, `tau = 2e-6`):

| `Delta` (per bp) | `Delta / tau` | Extra coverage ratio across 100 kb | Reads as |
|---|---|---|---|
| 0 | 0 | 1.00 | null |
| 1e-6 | 0.5 | 1.11 | barely detectable |
| 2e-6 | 1.0 | 1.22 | small |
| 3e-6 | 1.5 | 1.35 | moderate |
| 4.5e-6 | 2.25 | 1.57 | clear |
| 6e-6 | 3.0 | 1.82 | obvious (twice the reference rate) |

`Delta_std`, absolute `Delta`, `Delta / tau`, and the realized `Delta / tau_eff` are all
recorded per scenario, because the surface (C5) is probably indexed by the
last (2.5).

**Fixed or paired factors:**

| Axis | Main grid | Sensitivity values |
|---|---|---|
| Bin size | **5 kb** | 1 kb, paired on the same datasets (Q5) |
| Window | none (full genome) | 50 kb |
| Balance | balanced | 2 vs 5; 3 vs 8 (`G = 2`) |
| End anchors (S6) | none | anchor at 0 or at `L` (one arm) |

**Bin size.** 5 kb is the main grid and 1 kb the paired sensitivity: same
simulated atoms, one extra fit. The 1 kb package default is validated in full
for C1 and at the reference cells for C2. At `L = 50 kb`, 5 kb bins give only
10 bins per sample (and 1 or 2 bins on the short arm at anchor 0.1); the
short-genome cells are therefore expected to be fragile and are reported as
such, not dropped.

**Geometry and total decay.** Per-bp rates are held fixed across `L`, so a
shorter genome carries less total decay (`Delta * L`) and is harder; that is
the correct behavior for "my genome is smaller". Results are also tabulated
against the total log-decay across the longer arm so the two views reconcile.

### 4.2 Grid construction: a lean design justified by the collapse hypothesis

A full crossing of the design axes (7 x 3 x 6 = 126) with the geometry axes
(4 x 3 x 5 = 60) is 7,560 cells and is not run. Section 2.5 predicts that
`L`, `s_ref`, anchor, `sigma_atom`, and bin size act on power *only through*
`tau_eff`, so power collapses onto `(n, G, Delta / tau_eff)`. The grid is
built to test that prediction and then use it, rather than assume it:

- **Headline cells: reference geometry, full design crossing.** `L = 238,080`,
  `s_ref = 6e-6`, anchor 0.5, reference `tau` and `sigma_atom`; `n` x `G` x
  `Delta` in full (126 cells; Tier 1 subset 4 x 3 x 4 = 48). Headline Type I
  and power figures are reported here.
- **Geometry Type I map (Tier 1 and 2).** All 60 geometry combinations at
  `Delta = 0`, `G = 2`, `n` in {3, 5} (120 cells). Type I error does not
  depend on `tau_eff` in theory, so any geometry that departs from nominal is
  a finding and not explained by the collapse. This is the cheap, complete
  check on geometry.
- **Geometry power sweeps (one-at-a-time).** From the reference, vary one
  geometry axis at a time (3 + 2 + 4 = 9 variants) at `n` in {3, 5}, `G = 2`,
  `Delta` in {0, 3e-6, 6e-6} (54 cells). Plot power against `Delta` and
  against `Delta / tau_eff`; collapse is supported if the curves overlay in
  the second but not the first.
- **Collapse validation in the surface fit (Tier 3).** The surface design (a
  Sobol sequence) varies `n`, `G`, `Delta`, `tau`, `sigma_atom`, `L`, `s_ref`
  and anchor fraction jointly. The collapse is tested by whether a model in
  `(n, G, Delta / tau_eff)` alone has held-out error below the 0.03 target
  against the same model with the geometry covariates added. If geometry
  covariates improve held-out error materially, the calculator takes more
  inputs and the plan says so; the grid is not widened beforehand.
- **Note on 4.3.** The A2-2 prep prompt cites section 4.3; this plan has
  only 4.1 and 4.2, and the grid text was updated there.
- **Nuisance sweeps** for `cv` (formerly `tau`) and `sigma_atom` are one-at-a-time at the
  reference geometry, `n` in {3, 5}, `G = 2`.

Tiers:

- **Tier 1:** reference-geometry design crossing (48 cells), geometry Type I
  map (120 cells), geometry power sweeps (54 cells), nuisance sweeps (about
  14 cells), and each of S1 to S6 at `n` in {3, 5}, `G = 2`, `Delta` in {0,
  3e-6, 6e-6}. About 300 cells at R = 200.
- **Tier 2 (Type I):** reference geometry `Delta = 0`, `n` x `G` full (21
  cells) at 5 kb and 1 kb; partial nulls for `G` in {3, 4}; the geometry Type
  I map at higher R; S1 to S6 at reference cells chosen from Tier 1.
- **Tier 3 (power):** reference geometry full factorial (7 x 3 x 5 nonzero =
  105 cells) for the figures; about 400 Sobol points for the surface.

---

## 5. Infrastructure

### 5.1 Targets structure

The unit of work is a **batch**: one scenario times a fixed block of
replicates (default 100). Dynamic branching goes over batches, so branch count
is cells x R / 100 rather than cells x R.

```
scenario_spec (file target: scenarios by tier, tracked)   -> tar_target(format = "file")
grid           = build_grid(tier)                         # one row per scenario, with scenario_id and a hash
batches        = expand_batches(grid, R, batch_size)      # one row per (scenario_id, batch_id)
pkg_methods    = file target on the package's method files   # see below
sim_batch      = run_batch(batches, pkg_methods)          # pattern = map(batches); returns a tibble of per-replicate rows
sim_all        = bind_rows(sim_batch) per tier
summary_typeI  = summarize_typeI(sim_all)                 # rates, MC SE, Wilson/Clopper-Pearson CI, p-value ECDF
summary_power  = summarize_power(sim_all)
summary_bench  = summarize_paired(sim_all)                # paired LMM vs KW+Dunn
surface_fit    = fit_power_surface(summary_power)
figures        = make_figures(...)                        # file targets
tables         = write_tables(...)                        # file targets
```

Design choices worth flagging:

- **Seeds.** Each replicate gets its own seed derived from a hash of
  `(scenario_id, replicate_id)`, not a stream position. Any single failing
  replicate can be regenerated alone, extra replicates can be added later
  without disturbing existing ones, and results do not depend on how
  batches are scheduled across workers. Paired comparisons (baseline vs
  sensitivity, 1 kb vs 5 kb) reuse the same replicate seed for the shared
  random components. `null_distribution_F(seed = ...)` calls `set.seed`
  and would clobber the data-generating stream (parked item 10); the driver
  therefore generates the full dataset first and sets the seed explicitly
  per secondary-test call.
- **One simulation, several analyses.** A batch generates each dataset once
  and runs every analysis arm on it (LMM at 5 kb and 1 kb, KW + Dunn, the
  S2 drop-zero-bins variant, optional comparators), so arms are paired by
  construction and generation cost is not repeated.
- **Provenance targets.** `pkg_methods` is a file target covering only the
  files that define the statistics: `R/decay_tests.R`,
  `R/null_distribution.R`, `R/decay_slopes.R`, `R/coverage_io.R`, and the
  benchmark function once it exists. A change to any of these invalidates
  downstream simulation branches; unrelated package edits (plot code,
  vignettes, `man/`) do not trigger a rerun of the final grid. The package
  git SHA and `sessionInfo()` are recorded as a separate always-run target and
  written to the manifest, but are not a dependency of the simulation
  branches.
- **Failure handling.** `run_batch` wraps each analysis arm in `tryCatch` and
  returns a status column rather than dropping the row. Summaries report
  rates both unconditionally (failure counted as non-rejection) and
  conditional on success, plus the failure rate itself.
- **No nested parallelism.** Parallelism is at the branch level only. Inside a
  branch everything is serial (`BiocParallel::SerialParam()` for the
  secondary test), and BLAS/OpenMP threads are pinned to 1 so 10 workers do
  not oversubscribe a 12-core Mac.

### 5.2 Parallelism with `future`

As requested, the plan uses `future`: `future::plan(multisession, workers =
N)` set in `_targets.R`, run with `targets::tar_make_future()`. Two facts
from the installed `targets` 1.12.0 (checked with its help pages):
`tar_make_future()` is documented as **superseded** in favor of `tar_make()`
with the `crew` package, and `crew` and `crew.cluster` (which has a SLURM
controller, relevant for Discovery) are **not installed**, while `future` is.
The study builds with `future` now and keeps every target executor-agnostic
so switching to `crew` later is a change to `_targets.R` only. Adding `crew`
would be a new dependency and a separate decision; it is not in the `renv`
library (section 5.6).

### 5.3 Where intermediate results live

- The `targets` cache (`simulations/_targets/`) is local, large, and **not
  committed**. It is the memoization layer, not a deliverable.
- Committed outputs are written by file targets into `simulations/results/`
  and `simulations/figures/`.
- Per-replicate raw results for the *final* tiers only are committed as
  compressed `.rds` (base R, `xz`) so reviewers can recompute any summary
  without rerunning. Rough size: 0.5 to 1 million rows at about 100 bytes
  uncompressed each, so tens of MB before compression; to be measured in
  Tier 1. If it exceeds what is reasonable for git, the fallback is raw
  results on shared storage plus a checksum in the manifest. The simulation
  stack avoids `qs` and `arrow` (not installed) for the same dependency
  reason.
- `simulations/.gitignore` ignores `_targets/`, `renv/library/`,
  `renv/local/`, `renv/python/`, `renv/staging/`, and `results/scratch/`.
  `renv.lock` is committed. `.Rbuildignore` already excludes `simulations`.

### 5.4 File layout

```
simulations/
  PLAN.md                  this document
  README.md                how to rerun, tier by tier; written in the first execution session
  .gitignore               _targets/, renv/ library dirs, results/scratch/
  _targets.R               pipeline definition; selects tier via an environment variable
  renv.lock                pinned simulation dependencies (committed); renv/ activate script
  R/
    scenarios.R            grid definitions per tier, scenario_id hashing, seed derivation
    generate.R             generative model: atoms -> bins -> normalize_coverage -> tibble
    sensitivities.R        S1 to S5 as modifiers of the baseline generator
    analyze_lmm.R          wrapper around test_decay_slopes; extracts and stores 2.4 quantities
    analyze_bench.R        three benchmark arms (A2-2)
    analyze_secondary.R    wrapper around test_decay_significance (Tier 4)
    run_batch.R            replicate loop, error capture, timing
    summarize.R            rates, MC SE and CIs, paired differences, p-value ECDF
    surface.R              power surface fitting and inversion
    plot_theme.R           Okabe-Ito (categorical) and viridis (continuous) helpers
  scripts/
    pilot_timing.R         Tier 0
    check_generator.R      parameter recovery, atom-path vs per-base-path equivalence, NB validation
    run_local.R            tier launcher for the Mac (sets the future plan)
    run_hpc.sh             SLURM launcher, written only if Q4 resolves toward HPC
    make_figures.R, make_tables.R, fit_power_surface.R
  tests/
    testthat/              checks on the generator and summaries (see below)
  data/                    calibration inputs only (real-data summary statistics, Q1)
  results/
    explore/               Tier 1 summaries (committed)
    tables/                final tables
    surface/               fitted surface and its predictions
    raw/                   per-replicate results, final tiers only
    MANIFEST.csv           one row per committed artifact: tier, scenario family, R, package SHA, simulation SHA, date
  figures/                 final and supplementary figures
```

`tests/` (approved) holds the generator's own checks (see 6), kept out of the
package suite so package tests verify what the package claims and the
simulations verify that the claims are reasonable. They run with
`testthat::test_dir("simulations/tests")`, not `devtools::test()`.

### 5.5 Main driver

`_targets.R` reads one environment variable (`CDA_SIM_TIER`, values `pilot`,
`explore`, `typeI`, `power`, `secondary`), builds the grid for that tier, and
declares the pipeline in 5.1. `scripts/run_local.R` sets
`future::plan(multisession)`, pins BLAS threads, and calls `tar_make_future()`.
Re-running a tier after a partial failure resumes from the cache; raising `R`
for a scenario adds batches without recomputing existing ones.

### 5.6 Dependency pinning (`renv`)

Simulation-only dependencies are pinned by `simulations/renv.lock` (Q8), with
`renv` initialized scoped to `simulations/` (`renv::init(bare = TRUE)` run
from that directory), so the simulation library is separate from the package
and system libraries. The `.gitignore` in 5.3 keeps the library out of git and
the lockfile in.

---

## 6. Generator and summary validation (before any tier is trusted)

These are checks on the instrument, kept separate from the claims.

1. **Parameter recovery.** At large `n` (say 40 per group) and many bins, the
   fitted `Delta`, `tau`, and `sigma` recover the generating values.
2. **Idealized calibration.** With the baseline model and large `n`, Type I
   is near 0.05 and the p-value ECDF is near uniform. If it is not, the
   generator or the summary code is wrong before the method is.
3. **Positive control.** A very large `Delta` gives power near 1.
4. **Path equivalence.** Atom path and per-base `load_all_samples` path agree
   on `log_norm_cov`.
5. **Determinism.** Regenerating one replicate from its seed reproduces the
   stored result exactly; adding replicates leaves existing ones unchanged.
6. **Normalization is exercised.** Changing `D_i` alone does not change any
   test result.

---

## 7. Deliverables and naming

| Artifact | Path | Notes |
|---|---|---|
| Type I table (C1, C1b) | `results/tables/tbl-typeI-main.csv`, `tbl-typeI-posthoc.csv` | rate, MC SE, CI, n failed, singular rate, per scenario |
| Power table (C2) | `results/tables/tbl-power-main.csv`, `tbl-n-for-80pct.csv` | |
| Benchmark table (C3) | `results/tables/tbl-bench-lmm-vs-kwdunn.csv` | paired differences with CI |
| Sensitivity tables (C4) | `results/tables/tbl-sens-autocorr.csv`, `tbl-sens-dropout.csv`, `tbl-sens-hetero.csv`, `tbl-sens-anchor.csv` | change vs baseline |
| Secondary (C6) | `results/tables/tbl-secondary-typeI.csv`, `tbl-secondary-power.csv` | only if approved |
| Per-replicate raw | `results/raw/<tier>-<family>.rds` | final tiers only |
| Fitted power surface (C5) | `results/surface/power-surface-v1.rds` plus `power-surface-v1-grid.csv` (predictions on a regular grid) and `power-surface-v1-card.md` (domain of validity, held-out error, fit method, inputs the calculator needs) | versioned by suffix; a refit increments the version |
| Figures | `figures/fig-typeI-by-n.pdf`, `fig-typeI-pvalue-ecdf.pdf`, `fig-power-curves.pdf`, `fig-bench-power.pdf`, `fig-sens-*.pdf`, `fig-surface.pdf`; supplement as `figS-*.pdf` | PDF plus PNG; Okabe-Ito categorical, viridis continuous, no red-green |
| Pipeline | `simulations/_targets.R`, `R/`, `scripts/`, `renv.lock` | the reviewer-rerunnable artifact |
| Provenance | `results/MANIFEST.csv`, `results/sessionInfo.txt` | |

Figure numbers in file names are descriptive (`typeI-by-n`), not paper
numbers, so the Track B draft can renumber without renaming files.

**Handoff to the package.** The calculator function lives in the package, but
the surface artifact is produced here. A later package session copies
`power-surface-v1.rds` into `inst/extdata/` and writes the calculator against
its documented card. Nothing in this study writes to `inst/` or `R/`.

**Surface fitting method (proposal).** Two candidates compared by held-out
error: (a) a binomial GLM on the noncentral-F form, power = `1 - pf(qf(1 -
alpha, df1, df2), df1, df2, ncp)` with `ncp` modeled as a smooth function of
`n`, `G`, and standardized effect, which is monotone in effect size by
construction and interpretable; (b) a binomial GAM (`mgcv`, installed) on
`n`, `G`, standardized effect and nuisance covariates. (a) is preferred if it
fits within the held-out target, because the calculator must never return
power that decreases with more samples.

---

## 8. Proposed execution sessions

Mirroring the session structure in the package `PLAN.md`. Each ends with a
session log; none commits.

1. **A2-1: generator, validation, pilot.** `generate.R`, `scenarios.R`,
   `check_generator.R`, `analyze_lmm.R`, `pilot_timing.R`, `renv` and
   `_targets.R` skeleton, Tier 0. Output: measured time per fit, singular and
   failure rates, calls on M1 to M7, and a revised compute table. No arm
   beyond the primary LMM is implemented.
2. **A2-2: analysis arms and Tier 1.** `analyze_*.R`, `run_batch.R`, targets
   pipeline, Tier 1, collapse test (2.5), grid freeze proposal.
3. **A2-3: Tier 2 (Type I) and sensitivities S1 to S5.**
4. **A2-4: Tier 3 (power), benchmark, surface.**
5. **A2-5 (conditional): Tier 4 secondary test.**
6. **A2-6: figures, tables, manifest, README.**

Q1 to Q10 are resolved (section 9), so A2-1 is unblocked.

---

## 9. Decisions and risks

### Resolved (Carly's review, A2-1 start)

| Q | Decision |
|---|---|
| Q1 Calibration data | Superseded by scope: the study is not calibrated to HCMV. Parameters are centers of plausible ranges and the ranges are shown. No `simulations/data/` calibration inputs are required. |
| Q2 Benchmark spec | Holm for Dunn. Compute gated (KW `p < 0.05`) and ungated. Three variants: KW-chi2 + Dunn-Holm, KW-exact + Dunn-Holm, two-stage parametric (per-sample slope, then Welch's t for `G = 2` or one-way ANOVA for `G >= 3`). The exact variant is required (M4). Implemented in A2-2. |
| Q3 Secondary test (C6) | In scope at reduced size, Tier 4. Not implemented in A2-1. |
| Q4 Compute | Not decided; remains parked. The Tier 0 timing decides whether `run_hpc.sh` is written. |
| Q5 Bin size | 5 kb main grid; 1 kb paired sensitivity on the same datasets. |
| Q6 Failure counting | Headline is unconditional (failure counted as non-rejection); conditional-on-success reported alongside; the failure rate itself is stored. |
| Q7 Acceptance | Strict Clopper-Pearson. Type I is "controlled" in a cell if the upper end of its 95% CI at `alpha = 0.05` is at most 0.06; "liberal" if the lower end exceeds 0.05. Both recorded. |
| Q8 Pinning | `renv` lockfile at `simulations/renv.lock`, scoped to `simulations/`. |
| Q9 `G >= 3` pattern | One group differs from the shared reference (2.2). Not equal spacing. |
| Q10 Alpha | 0.05 for surface fitting. P-values are stored so other alphas are free in summaries. |

Approved additions: the three-arm benchmark, a random-slope-only comparator
conditional on Tier 0 confirming M1, and `simulations/tests/`.

**Risk: the generator is the study.** Results are statements about behavior
over the stated parameter ranges, with the ranges shown, not about any one
organism. The ranges for `tau`, `sigma_atom`, dropout rates and run lengths
are judgment-based and are the study's main external-validity assumption.

**Risk: dependency on package state.** The benchmark function does not exist
yet. Tier 2 and 3 results for C3 cannot be produced until it does, and any
later edit to the five method files (5.1) forces a rerun of the affected
tiers.

**Risk: failure modes dominate small `n`.** At `n = 2` per group the model has
4 samples and a 2 x 2 random-effects covariance, which M1 predicts is
degenerate. Expect a large share of the grid's low-`n` cells to be reported
with caveats rather than clean numbers.

---

## 10. Methodological observations from re-reading the protected files

Format per your instruction: I think X because Y, and Z would look like. None
of these changes protected files; the proposals are either simulation-arm
additions or decisions for you. M1 and M2 change the simulation design.

**M1. The random-effects covariance is probably degenerate by construction,
so `singular` may be `TRUE` for most fits and carry no information.
[to verify in Tier 0] Post-A2-1 decision, see `NOTES.md`.**
I think this because `normalize_coverage` divides each sample by its own
median bin coverage, so each sample's `log_norm_cov` is about 0 at the median
distance, whatever its depth. For a log-linear truth, the sample intercept at
the anchor is then determined by its slope (`c_i = s_i * d_median`), so
`b0_i` and `b1_i` in `(1 + distance | sample)` are almost perfectly
(negatively) correlated, which is the boundary of the unstructured 2 x 2
covariance; `isSingular` should fire. Shifting the distance origin does not
cure it: with an unstructured covariance the model is invariant to that shift,
and at the median distance the intercept variance goes to zero instead.
Consequences: (a) `print.decay_test` would warn "results should be
interpreted with caution" on nearly every fit, so the warning loses its
signal; (b) Kenward-Roger adjustments computed at a boundary estimate are the
case where their small-sample guarantee is weakest; (c) the Type I claim in C1
has to be reported conditional on, and unconditional of, `singular`.
What I would do: Tier 0 measures the singular rate and the correlation
estimate. If the singular rate is high, Tier 1 adds one comparator arm in the
simulation driver, `(0 + distance_scaled | sample)` (a random-slope-only model, no random
intercept), and compares Type I and power
with the package's model. If the comparator is better calibrated, that is a
methods memo to you with the evidence, not a change I would make.

**M2. The primary Type I claim depends on the between-sample slope variance,
and `n = 2` per group has few degrees of freedom for it.**
I think this because the interaction is tested against variation of
sample-level slopes, so the denominator df is at most about `N - G` (Kenward-
Roger will report it, stored as `df2`). At `G = 2`, `n = 2` that is 2 df;
`qf(0.95, 1, 2)` is 18.5, so the squared standardized slope difference has to
be large before the test rejects, however many bins each sample has. `R/decay_tests.R`'s
documentation describes `n = 2` to `5` as the target regime; the simulation
should say plainly what power exists there. Z: C2 reports power by `n` with
`df2` shown, and the surface card states a minimum `n` below which the
calculator refuses to return a number instead of returning a low one. This
also supports the 2.5 prediction that power is governed by `n` and
`Delta / tau`, not by bin count.

**M3. The anchor-specificity null is not exchangeable in the way the vignette
assumes, because F at an anchor depends on the anchor's geometry.**
I think this because `null_distribution_F` draws anchors uniformly over the
genome, but at an anchor near an end the fit uses one arm and a different
distance range (and `Sxx`) than at an interior anchor. The F reference
distribution under "nothing is special" therefore varies by position, so the
observed anchor (for example OriLyt, interior) is compared against a null
that mixes interior and end anchors. The vignette's exchangeability
assumption is stated for biological features, not for this geometric one.
A second contributor is redundancy: at 5 kb bins, anchors inside the same bin
produce nearly the same F, so `n_random = 500` on 238 kb is about 48 distinct
null values repeated about ten times each, and the `(k+1)/(n+1)` granularity
looks finer than it is (the roxygen already notes the reduced effective size).
Z: C6 simulates the specificity null with observed anchors at several
positions (interior and ends) and reports the empirical-p Type I rate by
anchor position, and the distinct-F count of the null. If the rate varies by
position, the alternative to propose is a null restricted to anchors of
matched geometry (for example, anchors with the same number of arms and a
similar distance range), written up as a memo, since that is a change to the
test machinery.

**M4. Kruskal-Wallis on per-sample slopes cannot reject at `n = 2`, and at `n
= 3` it rejects at about twice the nominal rate, because of the chi-square
approximation. [hand-computed, verify in Tier 0]**
I think this because with continuous slopes and `G = 2`, the KW statistic
takes only a few values. At `n = 2` per group the maximum (perfect separation)
is `H = 2.4`, chi-square(1) p = 0.12, so KW has zero power and zero Type I at
`n = 2` for `G = 2`; for `G = 3`, `n = 2`, the maximum `H` is 4.57, p = 0.10.
At `n = 3`, `G = 2`, only perfect separation (probability `2 / C(6,3) = 0.10`
under the null) gives `H = 3.857`, p = 0.0495, so the chi-square-approximate
test rejects in 10% of null datasets. The exact permutation test would give p
= 0.10 there and never reject. Consequently, the head-to-head at small `n` is
driven by discreteness and approximation, not by statistical efficiency.
Z: C3 reports KW both with the chi-square approximation and with exact
permutation p-values, and includes a third comparator arm, a two-stage
parametric analysis (per-sample slope, then a one-way ANOVA or Welch test),
so the LMM's advantage can be split into "mixed model beats two-stage" and
"parametric beats rank-based". The third arm is an addition for you to
approve; it is cheap.

**M5. A zero-coverage bin is the highest-leverage point in the model, and the
size of the problem depends on bin size, not just dropout rate.**
I think this because `log(0 + 1e-5)` is about -11.5 against typical
`log_norm_cov` near 0 to -2, and `test_decay_slopes` removes only non-finite
values. A single zero bin at the far end of a sample moves that sample's
slope estimate by an amount that dominates the between-sample slope variance
at `n = 2` to 3, and the random slope can absorb part of it, which inflates
`tau` and costs power, or fail to absorb it, which biases `Delta`. Because
binning averages atoms, a 5 kb bin reaches exactly zero far less often than a
1 kb bin at the same dropout rate, so bin size is partly a dropout-robustness
choice. Z: S2 measures this directly with an arm that drops `norm_cov == 0`
bins before the test. If dropping wins, the recommendation to you would be
documentation (a user-side filter in the intro vignette) or a package option,
and either is outside this study. I would not change the `1e-5` offset
without that evidence.

**M6. The `emmeans` steps run on every fit and dominate the call cost for the
simulation.**
I think this because `test_decay_slopes` always computes `emtrends` and
Tukey `pairs`, not only when the interaction is significant (the
documentation says "when the overall interaction is significant"). For
`null_distribution_F`, which uses only `F_stat`, this is wasted time on
every one of the `n_random` fits. For this study it is a constant-factor
cost, not a validity issue, and I have budgeted for it. If Tier 0 shows
`emmeans` takes more than a third of the call, the proposal would be an
internal helper that fits the model and returns the KR table only, with
`test_decay_slopes` and `null_distribution_F` both built on it. That is a
refactor of a protected file, so it is a request for permission, not
something I would do.

**M7. Smaller design notes from the code.**
- `test_decay_slopes` sets `check.conv.*` to `"ignore"`, so non-convergence
  is silent unless `conv_code` is read. The driver stores it, and summaries
  report Type I by convergence status.
- `test_decay_significance` with `window = NULL` defaults to no exclusion, so
  the null pool includes anchors a few kb from the observed one, which share
  its signal when the effect is localized. C6 power is therefore evaluated
  under both `window = NULL` and `window = 50 kb` (exclusion of 50 kb) to show
  the cost.
- The `n_samples < 3` and per-group `< 2` errors mean the smallest valid
  design is `G = 2`, `n = 2`; the grid starts there and no cell triggers
  those errors by construction.
