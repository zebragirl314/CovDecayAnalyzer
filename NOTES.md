# Project notes and decisions

Running record of scope decisions, deferred work, and methodology choices
for CovDecayAnalyzer. This file is intentionally terse. The vignettes are
the long-form documentation.

## Scope

CovDecayAnalyzer provides tools for analyzing coverage decay around a
user-supplied anchor position in linear reference coordinates. The method
is applicable to dsDNA viruses and other organisms with a localized
replication origin, whether the physical genome is linear or circular,
provided the reference sequence is represented linearly. Primary user
flow: read per-base coverage from `breseq`-style files or any equivalent
schema, bin and median-normalize per sample, fit log-linear decay slopes
from an anchor, test for group differences with a linear mixed-effects
model (Kenward-Roger F), and assess anchor-specificity against a random-
anchor null distribution.

Explicitly out of scope: reads-per-million quantification (lives in the
automated SLURM pipeline, not this package); per-feature coverage
differential analysis (use csaw or windowed DESeq2); coordinate wrap-
around (a user with a circular physical genome represented in a linear
reference can analyze their data today, but the method does not treat
position `genome_length` as adjacent to position `1`; this limitation is
addressed in the paper's Discussion).

## Methodological commitments

Primary test is the linear mixed-effects model with Kenward-Roger
small-sample correction, applied to all per-bin observations across all
samples. Rationale and formal specification are in
`vignettes/null-distribution.Rmd`.

Anchor-specificity test is a random-anchor Monte Carlo null distribution of
interaction F-statistics. The `(k+1)/(n+1)` empirical p-value is used as a
conservative convention; see the null-distribution vignette for the exact
framing.

The lab's older Kruskal-Wallis plus Dunn two-stage approach is superseded
and reimplemented in the simulation study as benchmark variants: KW with
chi-square p-values (traditional), KW with permutation p-values (correct
at small n), two-stage Welch, and two-stage pooled-variance ANOVA.

Carly Bobak holds editorial authority on method decisions for this package.

## Paper's methodological argument (post-A2-4)

The simulation study has tested the pre-simulation methodological argument
and refuted parts of it. The paper reports the simulation findings as the
scientific contribution rather than defending the original pivot narrative.

Honest statement of findings from Tier 2:

The primary LMM is well-calibrated across the parameter space tested.
Type I error ranges 0.043 to 0.054 across all 21 null cells of the main
grid at R=5000, including G=2 n=2 (0.045). Post-hoc contrasts control
FWER at or below 0.03 at partial-null cells.

The primary LMM and the pooled-variance two-stage ANOVA are
indistinguishable across balanced designs. In Stage 2 the two methods
agreed on all 63 figure-grid cells within 0.006 in paired rejection rate.
In Stage 3 S0 (unbalanced) the largest paired difference was 0.002. The
LMM does not have a general performance advantage over the pooled two-
stage test in the parameter space tested.

The LMM's small-n advantage at G=2 that was documented in Tier 1 is
specifically against Welch's t-test, not against the pooled two-stage
test. The pooled two-stage test ties the LMM at G=2 as well.

Under dropout rates of 2 to 5 percent at 5 kb run lengths, both the LMM
and the pooled two-stage are ultra-conservative. The LMM is more
conservative (Type I 0.005 to 0.011 at affected cells; pooled 0.020 to
0.043). The LMM's power loss relative to pooled is up to 0.112 at paired
comparisons on effect cells. Neither method is calibrated under these
conditions; the paper recommends preprocessing to address dropouts
before applying either test.

The slope-only parameterization (random slope with no random intercept)
has Type I of 0.507 (95% CI 0.493 to 0.521) at R=5000 reference and
shows bin-size-dependent inflation (+0.037 from 5 kb to 1 kb at R=1000).
The near-degenerate intercept-slope covariance in the full model is
absorbing real correlation structure that the Kenward-Roger correction
handles correctly.

The random-anchor null test (`test_decay_significance`) is a capability
the two-stage framework does not provide and is the primary unique
contribution of the package's statistical machinery. Simulation
validation of this test is in `vignettes/null-distribution.Rmd` and the
`R/null_distribution.R` test suite (coverage 99%).

What the simulation study did not establish:

Heteroscedastic residual variance (S3) as implemented could not
distinguish the three methods because the imposed noise term did not
reach the test at the reference sigma_atom. A targeted rerun at larger
sigma_atom is proposed as a follow-up before the paper's methods section
is finalized.

Heteroscedastic between-sample slope variance (S4) inflates Type I in
all three methods approximately equally (LMM, pooled, Welch). This is
not a methodological differentiator but is a documented condition where
all three methods misbehave similarly.

The approximate collapse hypothesis: power depends approximately on
`Delta/tau_eff` across s_ref levels, with deviations at the finest
resolution (+0.054 and +0.035 at G=2 Delta_std=0.5 n in {5, 10}) matching
the noncentral-F prediction with the explicit tau_eff that combines
between-sample slope variance and within-sample residual noise. The
sample size calculator indexes on `Delta/tau_eff` with that explicit
definition, not on `Delta_std/cv`.

## Simulation study

Lives in `simulations/` at the repo root, excluded from the package build
via `.Rbuildignore`. The simulation study is a scientific instrument on
par with the test suite, not a sandbox. Outputs are committed to version
control.

Simulations serve three ends:
1. Characterize the primary LMM's calibration and power under a parameter
   space spanning plausible scenarios for dsDNA viruses with localized
   replication origins. The study is not HCMV-calibrated; HCMV is one
   application point within the space.
2. Compare the primary LMM against the three KW+Dunn variants and the
   two two-stage parametric variants on the same simulated data.
3. Produce the empirical power surface that will fit the sample size
   calculator added as a new package function.

Simulation-design requirements (established across A2-1 through A2-4):

- Baseline generative model: dsDNA viral coverage structure with
  configurable anchor, effect size, n per group, bin size, additive
  noise, parameterized by `s_ref` as a scale parameter with between-
  sample variability as a coefficient of variation `cv = tau/s_ref` and
  group effect as a fractional change `Delta_std = Delta/s_ref`. The
  s_ref reparameterization was adopted at A2-2 after the original
  absolute-units parameterization confounded the collapse hypothesis.
- Group-effect pattern for G >= 3 is "one group differs from the shared
  reference." NOT equal spacing. Specifically: for G=3 with effect size
  Delta_std, two groups have decay rate s_ref and one has
  s_ref * (1 + Delta_std); analogously for G=4.
- Simulate before normalization: the generator draws coverage on a 1 kb
  atom grid, aggregates to the analysis bin size, and calls the
  package's own normalize_coverage(). Skipping normalization would
  change the random-effects structure and invalidate Type I estimates.
- Partial-matching convention for simulation code: all simulation-side
  function signatures should prevent R's partial argument matching from
  silently rewriting inputs. The A2-2 prep found `args$Delta` partial-
  matching `Delta_std`; use `match.arg()` or explicit validation.
- Integer overflow convention: any count-based computation that could
  exceed 2^31 at the full-R replicate counts should use numeric
  arithmetic by construction. The paired_diff() function hit this at
  A2-3 R=5000.

## Simulation findings parked for the paper

Headline numbers for the paper (all at reference nuisance unless noted):

- Primary LMM Type I curve: 0.043 to 0.054 across all 21 null cells of
  the main grid at R=5000. Full curve in `tier2/typeI/tbl-typeI-main.csv`.
- Primary LMM power at G=2, Delta_std=1, reference: 0.381 at n=2, 0.769
  at n=3, 0.985 at n=5. 95% CP CIs in `tier2/power/tbl-power-main.csv`.
- Smallest n for 80% power (LMM, G=2, reference): n=10 at Delta_std=0.5,
  n=4 at 1.0, n=3 at 1.5.
- 1 kb vs 5 kb at reference cell at R=1000: Type I 0.051 vs 0.049; paired
  power differences bounded within 0.009 across effect sizes.
- Anchor-position Type I at interior anchors across the geometric range
  tested: 0.0500 (0.044 to 0.056) at 0.25 L at R=5000. No inflation.
- Anchor at the extremes (0.0 and 1.0, single-arm geometry): soft-flag
  Type I straddling (0.063 to 0.070) at G=3 n=3. Mild but likely real.
- Dropout conservative behavior: Type I 0.005 to 0.011 at affected
  cells with 2 to 5 percent dropout rate at 5 kb runs. Power loss up to
  0.112 vs pooled two-stage at effect cells. Both methods below nominal
  but LMM more so.
- Permutation Dunn vs normal-approximation Dunn at G=3 n=5: agree within
  noise (0.962 vs 0.958). Dunn's chi-square approximation problem is
  specifically an omnibus problem, not a post-hoc problem.
- KW permutation Type I is zero at G=2 n<=3 by construction (smallest
  possible p-value is 0.10, which exceeds 0.05 for any n and G where the
  enumeration has fewer than 20 distinct values).

## Open design questions (post-A2-4)

- S3 heteroscedastic-residual rerun at sigma_atom in {0.3, 0.5}: the
  current S3 cannot distinguish the methods because the noise term does
  not reach the test. A 24-cell targeted rerun at larger sigma_atom
  would resolve whether the LMM has a true advantage under heteroscedastic
  residuals or ties the pooled two-stage there too. Decision pending.
- Dropout mechanism: the LMM's residual SD rises from 0.07 to 1.0 to 1.6
  under dropouts, the random-slope SD rises 1.5 to 3 times, and the
  singular-fit rate rises from 0.01-0.2 to 0.55-0.84. Why pooled two-
  stage is less affected is not yet tested. Investigation parked.
- Sample size calculator design: the collapse hypothesis is approximate
  under the full `Delta/tau_eff` definition; the surface fit session will
  decide whether to fit on `Delta/tau_eff` directly (cleaner derivation)
  or to fit on raw inputs with the collapse as a derived property
  (easier for users).

## Session summary convention

Every Claude Code session ends with a summary at
`session-logs/session-YYYY-MM-DD-<topic>.md` describing what changed, why,
any deviations from spec, open questions, and the state of
`devtools::test()` and `devtools::check()`. The summary is written for
review on a train or in bed, not for someone sitting in front of a diff.
It is staged with the code for Carly's review. See `PLAN.md` for the full
working-rules block.

## Punch list

Package work (Track A, substantially complete):
- Reclassify six exported helpers (DONE Session 4).
- Raise `R/null_distribution.R` test coverage above 80% (DONE Session 2).
- Bring vignettes up to date with `plot_sample` (DONE Session 3).
- Switch plot defaults to Okabe-Ito for colorblind-safe palettes
  (DONE Session 3).
- Promote the two-stage pooled-variance ANOVA and KW+Dunn variants from
  the simulation study to the package proper. Priority raised: the
  pooled two-stage is now the paper's "equivalent benchmark" and users
  should be able to call it from the package directly.
- Add sample size calculator function, fitted from simulation output.
  Collapse-approximate indexing on `Delta/tau_eff` with explicit noise-
  structure definition.
- `DESCRIPTION` text corrections (stale ANOVA, KW, permutation claims;
  stale RPM reference). Pending.
- `pbkrtest` must move from Suggests to Imports. Confirmed required at
  install time by the A2-1 clean renv install; `test_decay_slopes()`
  needs it for Kenward-Roger. Not blocking for A2 but blocking for any
  release.
- Reclassification follow-up pattern: future API-reclassification
  sessions should grep for usages of each demoted function in `inst/`,
  `vignettes/`, `README.md`, and `NEWS.md` and update them. Session 4
  missed the extdata README example and it was caught at A2-1.

Simulation study (Track A2):
- A2-1 through A2-4 complete. Tier 2 production run across Type I,
  power, and sensitivities finished. Next optional: S3 rerun at larger
  sigma_atom; required: power surface fit for the sample size calculator.

Methods paper (Track B):
- Methods document superseding the old KW doc on REDS: write after the
  Tier 2 commit lands. The methodological framing is now settled: LMM
  and pooled two-stage are equivalent across balanced and unbalanced
  designs; the LMM contributes the random-anchor null test as a unique
  capability; both methods are conservative under dropouts; the
  slope-only parameterization is a cautionary comparator.
- Methods paper draft follows the methods document.

Pipeline integration (Track C):
- Optional coverage-decay module in the SLURM pipeline, wrapping the
  package.
- Intake form field to toggle it.

Deferred:
- GitHub organization move. Package currently at
  `github.com/zebragirl314/CovDecayAnalyzer` (intern's personal account);
  destination TBD.
- Circular-biology-analyzed-as-linear sensitivity cell. Noted as a
  limitation in the paper's Discussion rather than tested.
- `plot_sample` fit line currently uses all sample bins, not only in-
  window bins. Revisit on next editorial pass of intro vignette.
- Simulation results size management: Tier 2 raw outputs at xz
  compression need measurement before Tier 2 commits.

## Baseline state

Tagged as `baseline-pre-claude` on 2026-10-22. At that tag: R CMD check
clean, all tests pass, 67% overall line coverage, docs regenerated with
roxygen2 8.1.1. This tag is an immovable reference point; use
`git diff baseline-pre-claude` to see cumulative changes from here.

Current state (uncommitted, post-A2-4): Track A complete (Sessions 1
through 4 committed), 158 tests passing, R CMD check clean with zero
notes, overall coverage above 91%. Track A2 Tier 1, A2-3, and Tier 2
Stages 1 through 3 complete and awaiting commit after the paper
framing is settled. Next: optional S3 rerun, then power surface fit.
