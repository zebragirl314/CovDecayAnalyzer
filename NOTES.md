# Project notes and decisions

Running record of scope decisions, deferred work, and methodology choices
for CovDecayAnalyzer. This file is intentionally terse. The vignettes are
the long-form documentation.

## Scope

CovDecayAnalyzer provides tools for analyzing coverage decay around a
user-supplied anchor position on a linear reference genome. Primary user
flow: read per-base coverage from `breseq`-style files or any equivalent
schema, bin and median-normalize per sample, fit log-linear decay slopes
from an anchor, test for group differences with a linear mixed-effects
model (Kenward-Roger F), and assess anchor-specificity against a
random-anchor null distribution.

Explicitly out of scope: reads-per-million quantification (lives in the
automated SLURM pipeline, not this package); per-feature coverage
differential analysis (use csaw or windowed DESeq2); circular genomes.

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
can compare the two routes head to head.

Carly Bobak holds editorial authority on method decisions for this package.

## Simulation study

Lives in `simulations/` at the repo root, excluded from the package build
via `.Rbuildignore`. The simulation study is a scientific instrument on
par with the test suite, not a sandbox. Outputs are committed to version
control.

Simulations serve three ends:
1. Validate Type I error and power of the primary test under realistic
   HCMV-like coverage structure.
2. Benchmark LMM plus Kenward-Roger against Kruskal-Wallis plus Dunn on
   the same simulated data.
3. Produce the empirical power surface that will fit the sample size
   calculator added as a new package function.

Simulation-design requirements (parked from the methodological observations
review):

- Baseline generative model: HCMV-like coverage structure with configurable
  anchor, effect size, n per group, bin size, additive noise.
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

## Session summary convention

Every Claude Code session ends with a summary at
`session-logs/session-YYYY-MM-DD-<topic>.md` describing what changed, why,
any deviations from spec, open questions, and the state of
`devtools::test()` and `devtools::check()`. The summary is written for
review on a train or in bed, not for someone sitting in front of a diff.
It is staged with the code for Carly's review. See `PLAN.md` for the full
working-rules block.

## Punch list

Package work:
- Reclassify six exported helpers (`bin_coverage`, `load_sample_coverage`,
  `normalize_coverage`, `read_coverage_tab` as internal;
  `add_distance_from_anchor` and `empirical_pvalue` kept exported).
- Raise `R/null_distribution.R` test coverage above 80% (currently 33%).
- Bring vignettes up to date with `plot_sample`, added in commit 22f7c25.
- Add Kruskal-Wallis plus Dunn benchmark function; name and multiplicity
  adjustment to be fixed in the simulation-study context.
- Add sample size calculator function, fitted from simulation output.
- Switch plot defaults to Okabe-Ito for colorblind-safe palettes.

Methodological fixes authorized 2026-10-22 (full detail in `PLAN.md`
Session 2):
- Items 1, 2 (plain bugs in `null_distribution_F`): resolved in Session 1,
  commit 58c7e2b.
- Item 3 (Phipson-Smyth documentation inconsistency): roxygen fix in
  Session 2; vignette prose drafted by Claude and edited by Carly.
- Item 4 (vignette "cannot occur by construction" claim): drafted by
  Claude and edited by Carly.
- Items 5, 6, 7, 9 (statistical and design decisions): code changes
  authorized for Session 2.
- Item 8: parked as a simulation-design requirement; see above.
- Item 10: minor polish for a later session.

Documentation and writing:
- Carly writes an updated methods document reflecting the LMM approach,
  to supersede the older Kruskal-Wallis document on REDS.
- Methods paper draft (Track B) with intro positioning LMM over KW,
  methods derived from the updated document, results from simulations plus
  eventual HCMV data.

Pipeline integration (Track C):
- Optional coverage-decay module in the SLURM pipeline, wrapping the
  package.
- Intake form field to toggle it.

Deferred:
- GitHub organization move. Package currently at
  `github.com/zebragirl314/CovDecayAnalyzer` (intern's personal account);
  destination TBD (BoscoResearchGroup already referenced in `DESCRIPTION`).
- `DESCRIPTION` text corrections (stale ANOVA, KW, permutation claims;
  stale RPM reference; `pbkrtest` Imports/Suggests question). Needs
  explicit `DESCRIPTION`-edit authorization. Grouped with a later release-
  hygiene session.

## Baseline state

Tagged as `baseline-pre-claude` on 2026-10-22. At that tag: R CMD check
clean, all tests pass, 67% overall line coverage, docs regenerated with
roxygen2 8.1.1. This tag is an immovable reference point; use
`git diff baseline-pre-claude` to see cumulative changes from here.
