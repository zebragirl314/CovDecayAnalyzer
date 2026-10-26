# Project notes and decisions

Running record of scope decisions, deferred work, and methodology choices
for CovDecayAnalyzer. This file is intentionally terse. The vignettes are
the long-form documentation.

## Scope

CovDecayAnalyzer provides a general framework for coverage decay analysis
around a user-supplied anchor position in linear reference coordinates.
The method is applicable to dsDNA viruses and other organisms with a
localized replication origin, whether the physical genome is linear or
circular, provided the reference sequence is represented linearly.
Primary user flow: read per-base coverage from `breseq`-style files or
any equivalent schema, bin and median-normalize per sample, fit log-
linear decay slopes from an anchor, test for group differences with the
recommended statistical test for the study design, and assess anchor-
specificity against a random-anchor null distribution.

Explicitly out of scope: reads-per-million quantification (lives in the
automated SLURM pipeline, not this package); per-feature coverage
differential analysis (use csaw or windowed DESeq2); coordinate wrap-
around (a user with a circular physical genome represented in a linear
reference can analyze their data today, but the method does not treat
position `genome_length` as adjacent to position `1`; this limitation is
addressed in the paper's Discussion).

## Paper framing (settled 2026-10-24)

The paper is a general framework introduction, not a methods-comparison
paper. The scientific contribution is the full coverage decay analysis
framework: the decay-slope-from-anchor framing, median-normalized binning
preprocessing, the random-anchor null for assessing anchor specificity,
the simulation infrastructure for method characterization, and the
recommendations for which test to use under which study design.

Audience: methods audience that reaches biologists working on dsDNA
viruses with localized replication origins. Primary example space is
HCMV, included as a case study using the lab's data; broader example
space is referenced in intro and discussion.

Position relative to prior work. The lab's previous ad-hoc application
of Kruskal-Wallis plus Dunn on per-sample slopes (used in a single
published paper, briefly mentioned in supplemental methods, never
formally written up) is a motivating case rather than a benchmark to
beat. The paper formalizes the framework that case was reaching for,
tests multiple statistical options against each other and against the
primary LMM, and provides evidence-based recommendations.

The three comparator roles in the paper:
- Primary LMM (`test_decay_slopes`): the package's recommended single-
  step test. Characterized across the parameter space.
- Pooled-variance two-stage ANOVA on per-sample slopes: parametric
  secondary comparator. Simulations show it is approximately
  equivalent to the primary LMM across the parameter space tested.
  Added to the package as a comparator function.
- Kruskal-Wallis plus Dunn on per-sample slopes: the method used
  previously by the lab, in three variants (chi-square p-values,
  permutation p-values, with the chi-square variant miscalibrated at
  small n and the permutation variant having a power floor at the
  smallest n). Documented as a historically-used alternative, not
  recommended for small n.
- Slope-only LMM (random slope, no random intercept): cautionary
  comparator demonstrating the random intercept is doing real work
  despite the near-degenerate covariance estimate.

Note on normality assumptions. The two-stage comparator tests per-sample
slope estimates (OLS coefficients over 20 to 50 binned observations),
not raw coverage values. By the usual central limit argument, slope
estimates are approximately normal under the generative model even when
coverage is not. Simulation validates this empirically: pooled two-stage
is calibrated across the parameter space tested at n as small as 2 per
group. The paper is explicit that the pooled two-stage operates on slope
estimates, not on raw coverage.

The random-anchor null test (`test_decay_significance`) is a capability
unique to this package in the comparator set; the two-stage framework
and KW+Dunn do not provide an equivalent. This is a primary novel
contribution.

## Recommendation table (post-Tier 2 findings)

Organized by study design, which is how users will consult the paper's
recommendation table.

Balanced designs, n per group >= 3, no known heteroscedasticity, low
dropout rate: primary LMM (`test_decay_slopes`). Pooled two-stage ANOVA
is an approximately equivalent alternative.

Unbalanced designs: primary LMM or pooled two-stage ANOVA (simulations
show they agree within 0.002 in paired rejection rate on all unbalanced
cells tested).

Known or suspected heteroscedastic residual variance between groups, at
small n (n=3): pooled two-stage ANOVA. The primary LMM is calibrated but
conservative, with up to 0.041 paired rejection difference below pooled
at Tier 2 S3 rerun cells. At n>=5 the two methods agree within 0.013.

High dropout rates (>= 2% at 5 kb run lengths): preprocess to address
dropouts before applying either test. Both the primary LMM and the
pooled two-stage are below nominal in Tier 2 (LMM Type I 0.005-0.011;
pooled 0.020-0.043). The LMM power loss relative to pooled is up to
0.112 at effect cells. Neither method is calibrated under these
conditions.

Smallest-n designs (G=2, n=2): primary LMM. Welch's two-stage t loses
power substantially (0.17 vs 0.39 for pooled at n=3 Delta_std=0.5); KW
permutation has a power floor (minimum p=0.10). Primary LMM achieves
0.381 power at n=2 Delta_std=1.

Anchor near the ends of the reference (0.0 L or 1.0 L, single-arm
geometry): mild Type I elevation (0.063-0.070 at G=3 n=3). All methods
affected similarly.

Assessing whether a signal at an anchor is specific to that anchor
versus a general feature of the reference: `test_decay_significance`
with the random-anchor null. Unique to this package.

## Methodological commitments

Primary test is the linear mixed-effects model with Kenward-Roger
small-sample correction, applied to all per-bin observations across all
samples. Rationale and formal specification are in
`vignettes/null-distribution.Rmd`.

Anchor-specificity test is a random-anchor Monte Carlo null distribution
of interaction F-statistics. The `(k+1)/(n+1)` empirical p-value is used
as a conservative convention; see the null-distribution vignette for the
exact framing.

Carly Bobak holds editorial authority on method decisions for this
package.

## Simulation study

Lives in `simulations/` at the repo root, excluded from the package
build via `.Rbuildignore`. The simulation study is a scientific
instrument on par with the test suite, not a sandbox. Outputs are
committed to version control.

Simulations serve three ends:
1. Characterize the primary LMM's calibration and power under a
   parameter space spanning plausible scenarios for dsDNA viruses with
   localized replication origins. The study is not HCMV-calibrated; HCMV
   is one application point within the space.
2. Compare the primary LMM against four comparators (slope-only
   cautionary, three KW+Dunn variants, two two-stage parametric
   variants) on the same simulated data.
3. Produce the empirical power surface that will fit the sample size
   calculator added as a new package function.

Simulation-design requirements (established across A2-1 through A2-4):

- Baseline generative model: dsDNA viral coverage structure with
  configurable anchor, effect size, n per group, bin size, additive
  noise, parameterized by `s_ref` as a scale parameter with between-
  sample variability as a coefficient of variation `cv = tau/s_ref` and
  group effect as a fractional change `Delta_std = Delta/s_ref`.
- Group-effect pattern for G >= 3 is "one group differs from the shared
  reference." NOT equal spacing.
- Simulate before normalization: the generator draws coverage on a 1 kb
  atom grid, aggregates to the analysis bin size, and calls the
  package's own normalize_coverage().
- Partial-matching convention: all simulation-side function signatures
  prevent R's partial argument matching from silently rewriting inputs.
- Integer overflow convention: count-based computation uses numeric
  arithmetic.

## Simulation findings parked for the paper

Headline numbers at reference nuisance unless noted:

- Primary LMM Type I curve: 0.043 to 0.054 across all 21 null cells of
  the main grid at R=5000. Full curve in
  `tier2/typeI/tbl-typeI-main.csv`.
- Primary LMM power at G=2, Delta_std=1, reference: 0.381 at n=2, 0.769
  at n=3, 0.985 at n=5. 95% CP CIs in `tier2/power/tbl-power-main.csv`.
- Smallest n for 80% power (LMM, G=2, reference): n=10 at Delta_std=0.5,
  n=4 at 1.0, n=3 at 1.5.
- 1 kb vs 5 kb at reference cell at R=1000: Type I 0.051 vs 0.049;
  paired power differences bounded within 0.009 across effect sizes.
- Anchor-position Type I at interior anchors: 0.0500 (0.044 to 0.056)
  at 0.25 L at R=5000. No inflation across the geometric range tested.
- Anchor at the extremes (0.0 and 1.0): soft-flag Type I straddling
  (0.063 to 0.070) at G=3 n=3. Mild but likely real.
- LMM vs pooled two-stage across balanced designs (Stage 2, 63 cells
  at R=1000): paired rejection differences within 0.006. Equivalent
  across the balanced parameter space.
- LMM vs pooled two-stage at unbalanced cells (Stage 3 S0, 30 cells):
  largest paired difference 0.002. Equivalent at unbalanced designs.
- LMM vs pooled two-stage under heteroscedastic residuals (S3 rerun at
  sigma_atom in {0.3, 0.5}, with sigma_ratio=2 vs baseline ratio=1):
  divergence specifically at n=3 with heteroscedasticity plus high
  noise. Largest paired difference 0.041. At n=5 the methods agree
  within 0.013. At homoscedastic residuals the methods agree within
  0.015.
- Dropout conservative behavior: Type I 0.005 to 0.011 at affected
  cells with 2 to 5 percent dropout rate at 5 kb runs. Power loss up
  to 0.112 vs pooled two-stage at effect cells. Both methods below
  nominal but LMM more so.
- Singular-fit mechanism: three different stressors (dropouts at 2-5%,
  heteroscedastic residuals ratio=2 at high sigma, high residual noise
  alone) all raise the singular-fit rate above 0.5 and all produce LMM
  conservativeness in the same direction. Diagnostic: when the random-
  effects covariance structure is stressed, the LMM widens uncertainty
  rather than make false claims. Feature, not bug.
- Slope-only cautionary: Type I 0.507 (CI 0.493-0.521) at R=5000
  reference; bin-size-dependent (+0.037 from 5 kb to 1 kb at R=1000).
- Permutation Dunn vs normal-approximation Dunn at G=3 n=5: agree
  within noise (0.962 vs 0.958). Dunn's chi-square approximation
  problem is specifically an omnibus problem, not a post-hoc problem.
- KW permutation Type I is zero at G=2 n<=3 by construction (smallest
  possible p-value is 0.10).
- Approximate collapse hypothesis: power depends approximately on
  `Delta/tau_eff` with deviations matching the noncentral-F prediction.
  Sample size calculator indexes on `Delta/tau_eff` with explicit
  tau_eff definition from the generative model.

## Open design questions (post-A2-4)

- Sample size calculator design: fit the power surface on `Delta/tau_eff`
  directly (cleaner derivation) or on raw inputs with the collapse as a
  derived property (easier for users). Decision at the surface fit
  session.
- Dropout mechanism investigation: why the pooled two-stage is less
  conservative than the LMM under dropouts. The fitted-value diagnostics
  suggest the OLS first stage downweights the dropout zeros naturally,
  but this has not been tested formally. Parked for the paper's
  discussion or follow-up work.

## Session summary convention

Every Claude Code session ends with a summary at
`session-logs/session-YYYY-MM-DD-<topic>.md`. The summary is written for
review on a train or in bed, not for someone sitting in front of a diff.
It is staged with the code for Carly's review.

## Punch list

Package work (Track A, substantially complete):
- Reclassify six exported helpers (DONE Session 4).
- Raise `R/null_distribution.R` test coverage above 80% (DONE Session 2).
- Bring vignettes up to date with `plot_sample` (DONE Session 3).
- Switch plot defaults to Okabe-Ito (DONE Session 3).
- Promote the pooled-variance two-stage ANOVA to the package. Priority
  raised: the paper recommends it as the equivalent benchmark and under
  heteroscedastic small-n designs. Users should be able to call it from
  the package directly.
- Promote the KW+Dunn variants to the package. Historically-used method,
  retained for reproducibility of prior lab analyses.
- Add sample size calculator function, fitted from simulation output.
  Indexed on `Delta/tau_eff` with explicit noise-structure definition
  from the generative model.
- `DESCRIPTION` text corrections (stale ANOVA, KW, permutation claims;
  stale RPM reference). Pending.
- `pbkrtest` must move from Suggests to Imports (confirmed required at
  install time by A2-1 clean renv install). Blocking for any release.
- Repo size management: Tier 2 commit landed at 70.58 MiB and future
  tier results will grow proportionally. Decision pending between
  xz-compressing per-replicate .rds, moving raw results to REDS with
  checksummed manifest, or Git LFS. Recommendation is REDS with
  manifest.

Simulation study (Track A2):
- A2-1 through A2-4 complete including S3 rerun and baseline. Next:
  power surface fit for sample size calculator.

Methods paper (Track B):
- Framing memo (methodological positioning, recommendation table
  structure, three-comparator distinction): Claude drafts, Carly edits.
- Methods document superseding Noelle's KW doc: Claude drafts, Carly
  edits. Framing is now settled.
- Methods paper draft: structure is settled (general framework paper,
  by-design recommendation table, HCMV as case study, simulation
  evidence throughout).

Pipeline integration (Track C):
- Optional coverage-decay module in the SLURM pipeline, wrapping the
  package.
- Intake form field to toggle it.

Deferred:
- GitHub organization move. Package currently at
  `github.com/zebragirl314/CovDecayAnalyzer`; destination TBD.
- Circular-biology-analyzed-as-linear sensitivity cell. Noted as
  limitation in paper Discussion rather than tested.
- `plot_sample` fit line uses all sample bins, not only in-window bins.
  Revisit on next editorial pass of intro vignette.

## Baseline state

Tagged as `baseline-pre-claude` on 2026-10-22. At that tag: R CMD check
clean, all tests pass, 67% overall line coverage, docs regenerated with
roxygen2 8.1.1.

Current state (committed as 55d8fcd on 2026-10-24): Track A complete
(Sessions 1 through 4 committed), 158 tests passing, R CMD check clean
with zero notes, overall coverage above 91%. Track A2 Tier 1, A2-3,
Tier 2 Stages 1 through 3, S3 rerun, S3 baseline rerun all committed.
Simulation study substantively complete. Next: NOTES.md commit, framing
memo, methods doc draft, power surface fit session.
