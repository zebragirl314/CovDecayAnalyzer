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
in the lab's older methods document, is superseded. It is reimplemented
in the simulation study as a benchmark function in three variants: KW
with chi-square p-values (traditional), KW with permutation p-values
(correct at small n), and a two-stage parametric alternative (per-sample
slope then Welch's or pooled t for G=2, one-way ANOVA for G>=3).

Carly Bobak holds editorial authority on method decisions for this package.

## Paper's methodological argument (post-A2-3)

The simulation study has shifted the paper's framing. The LMM's value is
specifically in three areas, not in a general "beats two-stage" claim:

- Unbalanced designs. Classical ANOVA on per-sample slopes assumes equal
  n per group; the LMM accommodates imbalance by construction.
- Heteroscedastic designs (residual variance or between-sample slope
  variance differing across groups). The LMM can be extended to model
  this; the two-stage framework cannot.
- The random-anchor null test, which the two-stage framework does not
  support.

At balanced designs with G>=3, the pooled-variance two-stage ANOVA ties
the primary LMM on every validation cell at Tier 1. The paper reports
this finding directly rather than overstating the LMM's advantage.

At G=2 with balanced designs, the LMM's small-n gain over the two-stage
is specifically about pooled vs Welch variance (not about the mixed-model
framework). The two-stage with pooled variance also ties the LMM at
G=2 n>=3.

The paper's cautionary paragraph will describe the slope-only
parameterization (random slope, no random intercept): Type I of 0.507
(95% CI 0.493-0.521) at R=5000 reference, confirming the random
intercept absorbs real correlation structure despite the near-degenerate
covariance estimate. The slope-only arm also shows bin-size-dependent
Type I inflation (+0.037 from 5 kb to 1 kb at R=1000), strengthening
the argument.

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
2. Benchmark LMM plus Kenward-Roger against three KW+Dunn variants and
   two two-stage parametric variants on the same simulated data.
3. Produce the empirical power surface that will fit the sample size
   calculator added as a new package function.

Simulation-design requirements (parked from the methodological
observations review):

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
- Autocorrelated residuals sensitivity (parked for Tier 2): residuals
  correlated across adjacent bins per sample.
- Zero-coverage dropout sensitivity (parked for Tier 2): dropouts
  generated as runs on the atom grid.
- Partial-matching convention for simulation code: all simulation-side
  function signatures should prevent R's partial argument matching from
  silently rewriting inputs. The A2-2 prep found `args$Delta` partial-
  matching `Delta_std`; use `match.arg()` or explicit validation.
- Integer overflow convention: any count-based computation that could
  exceed 2^31 at the full-R replicate counts should use numeric
  arithmetic by construction. The paired_diff() function hit this at
  A2-3 R=5000.

## Simulation findings parked for the paper's Discussion

- Collapse hypothesis supported at Tier 1: power at matched Delta_std/cv
  agrees across s_ref levels within Monte Carlo CI. The sample size
  calculator can be indexed by the standardized effect alone.
- Anchor-position Type I resolved at Tier 1 plus A2-3: no inflation at
  interior anchors across the geometric range tested. The 0.105 Tier 1
  reading at 0.25 L was a 200-replicate noise artifact; at R=5000 it is
  0.0500 (0.044-0.056).
- 1 kb vs 5 kb at reference cell at R=1000: Type I 0.051 vs 0.049,
  paired power differences bounded within 0.009 across effect sizes.
  No bin-size calibration issue.
- Permutation Dunn vs normal-approximation Dunn at G=3 n=5: agree
  within noise (0.962 vs 0.958). Dunn's chi-square approximation
  problem is specifically an omnibus problem, not a post-hoc problem.
- KW permutation Type I is zero at G=2 n<=3 by construction (smallest
  possible p-value is 0.10, which exceeds 0.05 for any n and G where
  the enumeration has fewer than 20 distinct values).

## Open design questions (post-Tier 2)

- M1: random-effects covariance degeneracy, resolved for now. 13%
  singular rate with intercept-slope correlation near -1 in most fits.
  Session 2 isSingular surface decision stands. The slope-only
  comparator (0.507 Type I) demonstrates the intercept is doing real
  work. No further investigation planned unless Tier 2 surfaces new
  evidence.
- M3: anchor at the extremes (0.0 and 1.0, single-arm geometry) not yet
  tested at R=5000. Carry into Tier 2 as a sensitivity cell.

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
- Add KW+Dunn benchmark function (three variants) to the package. The
  simulation study prototypes all three; promote to the package after
  Tier 2 confirms the final API. The pooled-variance two-stage is also
  a candidate for promotion.
- Add sample size calculator function, fitted from simulation output
  (A2-4 or later). Collapse hypothesis supports a one-input calculator
  indexed by standardized effect size.
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
- A2-1 through A2-3 complete. Tier 1 exploratory and three focused
  high-R cells committed. Next: A2-4 (Tier 2 proper).

Methods paper (Track B):
- Methods document superseding the old KW doc on REDS: write after
  A2-3 is committed. The two-stage parity finding reshaped the paper's
  argument; writing the methods doc now would require redoing it
  against the pivot.
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
- Simulation results size management: Tier 1 at 22,400 datasets was
  3 MiB; Tier 2 at ~400,000 datasets could reach 50 MiB. Options when
  needed: `xz` compression on per-replicate `.rds`, external storage
  with checksums committed, or Git LFS. Revisit before Tier 2 commits.

## Baseline state

Tagged as `baseline-pre-claude` on 2026-10-22. At that tag: R CMD check
clean, all tests pass, 67% overall line coverage, docs regenerated with
roxygen2 8.1.1. This tag is an immovable reference point; use
`git diff baseline-pre-claude` to see cumulative changes from here.

Current state (post-A2-3): Track A complete (Sessions 1 through 4
committed), 158 tests passing, R CMD check clean with zero notes,
overall coverage above 91%. Track A2 at Tier 1 exploratory plus three
focused high-R cells complete. Next: A2-4 (Tier 2 proper).
