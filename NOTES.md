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

The previously-used Kruskal-Wallis plus Dunn two-stage approach, described
in the lab's older methods document, is superseded. It will be re-
implemented in the package as a benchmark function so the simulation study
can compare the two routes head to head. Three variants of the benchmark
will be implemented: KW with chi-square p-values (traditional), KW with
exact permutation p-values (correct at small n), and a two-stage parametric
alternative (per-sample slope then Welch/ANOVA).

Carly Bobak holds editorial authority on method decisions for this package.

## Simulation study

Lives in `simulations/` at the repo root, excluded from the package build
via `.Rbuildignore`. The simulation study is a scientific instrument on
par with the test suite, not a sandbox. Outputs are committed to version
control.

Simulations serve three ends:
1. Validate Type I error and power of the primary test under a parameter
   space spanning plausible scenarios for dsDNA viruses with localized
   replication origins. The study is not HCMV-calibrated; HCMV is one
   application point within the space.
2. Benchmark LMM plus Kenward-Roger against KW plus Dunn (chi-square),
   KW plus Dunn (exact), and a two-stage parametric alternative on the
   same simulated data.
3. Produce the empirical power surface that will fit the sample size
   calculator added as a new package function.

Simulation-design requirements (parked from the methodological
observations review):

- Baseline generative model: dsDNA viral coverage structure with
  configurable anchor, effect size, n per group, bin size, additive
  noise, parameterized to span biologically plausible ranges for linear-
  reference-coordinate analysis of viral genomes.
- Autocorrelated residuals sensitivity. The LMM assumes residual
  independence given random intercept and slope. Adjacent bins from the
  same sample have correlated residuals when the true decay is not log-
  linear, which can make within-sample F anti-conservative. The direction
  of the effect on the between-sample interaction test is unclear and
  worth measuring directly.
- Zero-coverage dropout sensitivity. `normalize_coverage` adds 1e-5 before
  taking the log, so a zero-coverage bin produces a value near -11.5, a
  high-leverage outlier on the log scale. The simulation should test
  whether realistic dropout rates bias slope estimates.
- Group-effect pattern for G >= 3 is "one group differs from the shared
  reference." Equal spacing is not the default.

## Open design questions (post-Tier-0)

- M1: random-effects covariance degeneracy. Reading `test_decay_slopes`
  suggests that median normalization pins each sample's intercept to its
  slope, so the random-effects covariance may be degenerate by
  construction and `isSingular(fit)` may fire on nearly every fit. If
  Tier 0 confirms this prediction, the surface-isSingular decision from
  Session 2 needs revisiting. Three resolutions on the table: suppress
  the warning and leave `isSingular` as a programmatic field only; report
  the random-effects correlation estimate instead of just singularity; or
  switch the primary model to random-slope-only. All three are post-A2-1
  decisions.

## Session summary convention

Every Claude Code session ends with a summary at
`session-logs/session-YYYY-MM-DD-<topic>.md` describing what changed, why,
any deviations from spec, open questions, and the state of
`devtools::test()` and `devtools::check()`. The summary is written for
review on a train or in bed, not for someone sitting in front of a diff.
It is staged with the code for Carly's review. See `PLAN.md` for the full
working-rules block.

## Punch list

Package work (Track A):
- Reclassify six exported helpers (DONE Session 4).
- Raise `R/null_distribution.R` test coverage above 80% (DONE Session 2).
- Bring vignettes up to date with `plot_sample` (DONE Session 3).
- Switch plot defaults to Okabe-Ito for colorblind-safe palettes
  (DONE Session 3).
- Add Kruskal-Wallis plus Dunn benchmark function (three variants).
  Prototype in `simulations/` during A2-2; promote to package after
  simulation design settles the API.
- Add sample size calculator function, fitted from simulation output.
- `DESCRIPTION` text corrections (stale ANOVA, KW, permutation claims;
  stale RPM reference; `pbkrtest` Imports/Suggests question). Grouped
  with a later release-hygiene session.

Methodological fixes authorized 2026-10-22 (full detail in `PLAN.md`
Session 2): all applied in Sessions 1, 2.

Simulation study (Track A2):
- Open at A2-1 with generator, infrastructure, and Tier 0 pilot. See
  `simulations/PLAN.md` for execution sequence.

Methods paper (Track B):
- Carly writes an updated methods document reflecting the LMM approach,
  to supersede the older Kruskal-Wallis document on REDS. Position the
  method as general (dsDNA viruses with localized replication origins)
  with HCMV as a demonstration application.
- Methods paper draft with intro positioning LMM over KW, methods
  derived from the updated document, results from simulations plus HCMV
  demonstration section.

Pipeline integration (Track C):
- Optional coverage-decay module in the SLURM pipeline, wrapping the
  package.
- Intake form field to toggle it.

Deferred:
- GitHub organization move. Package currently at
  `github.com/zebragirl314/CovDecayAnalyzer` (intern's personal account);
  destination TBD (BoscoResearchGroup already referenced in `DESCRIPTION`).
- Circular-biology-analyzed-as-linear sensitivity cell. Noted as a
  limitation in the paper's Discussion rather than tested in the
  simulation study; most users analyzing HCMV already do exactly this
  and it works for them.
- Design-question `plot_sample` fit line uses all sample bins, not only
  in-window bins. Revisit on next editorial pass of intro vignette.

## Baseline state

Tagged as `baseline-pre-claude` on 2026-10-22. At that tag: R CMD check
clean, all tests pass, 67% overall line coverage, docs regenerated with
roxygen2 8.1.1. This tag is an immovable reference point; use
`git diff baseline-pre-claude` to see cumulative changes from here.

Current state: Track A complete (Sessions 1 through 4 committed), 158
tests passing, R CMD check clean with zero notes, overall coverage above
91%, Track A2 planning committed. Next: A2-1 execution.
