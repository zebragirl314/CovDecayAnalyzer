# Project notes and decisions

Running record of scope decisions, deferred work, and methodology choices for
CovDecayAnalyzer. This file is intentionally terse. The vignettes are the
long-form documentation.

## Scope

CovDecayAnalyzer provides tools for analyzing coverage decay around a
user-supplied anchor position on a linear reference genome. Primary user flow:
read per-base coverage from `breseq`-style files or any equivalent schema, bin
and median-normalize per sample, fit log-linear decay slopes from an anchor,
test for group differences with a linear mixed-effects model (Kenward-Roger F),
and assess anchor-specificity against a random-anchor null distribution.

Explicitly out of scope: reads-per-million quantification (lives in the
automated SLURM pipeline, not this package); per-feature coverage differential
analysis (use csaw or windowed DESeq2); circular genomes.

## Methodological commitments

Primary test is the linear mixed-effects model with Kenward-Roger small-sample
correction, applied to all per-bin observations across all samples. Rationale
and formal specification are in `vignettes/null-distribution.Rmd`.

Anchor-specificity test is a random-anchor Monte Carlo null distribution of
interaction F-statistics, with the Phipson-Smyth (k+1)/(n+1) empirical p-value.

The previously-used Kruskal-Wallis plus Dunn two-stage approach, described in
the lab's older methods document, is superseded. It will be re-implemented in
the package as a benchmark function so the simulation study can compare the
two routes head to head.

Carly Bobak holds editorial authority on method decisions for this package.

## Simulation study

Lives in `simulations/` at the repo root, excluded from the package build via
`.Rbuildignore`. The simulation study is a scientific instrument on par with
the test suite, not a sandbox. Outputs are committed to version control.

Simulations serve three ends:
1. Validate Type I error and power of the primary test under realistic
   HCMV-like coverage structure.
2. Benchmark LMM plus Kenward-Roger against Kruskal-Wallis plus Dunn on the
   same simulated data.
3. Produce the empirical power surface that will fit the sample size
   calculator added as a new package function.

## Punch list

Package work:
- Reclassify seven exported helpers (`add_distance_from_anchor`, `bin_coverage`,
  `empirical_pvalue`, `load_sample_coverage`, `normalize_coverage`,
  `read_coverage_tab`, and reassess as internal versus exported).
- Raise `R/null_distribution.R` test coverage from current 33% to above 80%.
- Bring vignettes up to date with `plot_sample`, the single-sample view added
  in commit 22f7c25.
- Add Kruskal-Wallis plus Dunn benchmark function. Likely exported for
  reproducibility; final name to be decided.
- Add sample size calculator function, fitted from simulation output. New
  user-facing API.

Documentation and writing:
- Carly writes an updated methods document reflecting the LMM approach, to
  supersede the older Kruskal-Wallis document on REDS.
- Methods paper draft (separate track) with intro positioning LMM over KW,
  methods derived from the updated document, results from simulations plus
  eventual HCMV data.

Pipeline integration (separate track):
- Optional coverage-decay module in the SLURM pipeline, wrapping the package.
- Intake form field to toggle it.

Deferred:
- GitHub organization move. Package currently at
  `github.com/zebragirl314/CovDecayAnalyzer` (intern's personal account);
  destination TBD (BoscoResearchGroup already referenced in `DESCRIPTION`).

## Baseline state

Tagged as `baseline-pre-claude` on 2026-10-22. At that tag: R CMD check
clean, all tests pass, 67% overall line coverage, docs regenerated with
roxygen2 8.1.1. This tag is an immovable reference point; use
`git diff baseline-pre-claude` to see cumulative changes from here.
