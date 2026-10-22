# PLAN.md: package work sequence for CovDecayAnalyzer

Status: proposal reviewed and approved for execution. Updated 2026-10-22
after Session 1 (bug fixes to `null_distribution_F`) landed and the
methodological observations were reviewed.

Scope: the "Package work" section of the `NOTES.md` punch list. The
simulation study (Track A2), the methods paper (Track B), and pipeline
integration (Track C) are excluded and will have their own planning
documents when opened.

## Working rules for every session

- Session-start routine from `CLAUDE.md`: read `NOTES.md`, `git log`,
  `git status`, confirm the task still fits before any edit.
- Working-tree rules: uncommitted changes on `main`, no commits, no
  branches, no pushes.
- Write permissions in effect for the sessions in this plan: `R/`,
  `tests/`, `vignettes/`, `simulations/`, `inst/`, `NOTES.md`, `NEWS.md`,
  `README.md`, `.Rbuildignore`, `session-logs/`. `DESCRIPTION`, `NAMESPACE`,
  and `CLAUDE.md` remain protected without explicit instruction.
- Protection on `R/decay_tests.R` and `R/null_distribution.R` is lifted
  only for the specific authorized changes documented in each session
  block below. Other edits to those files remain protected.
- Every session ends with a summary at
  `session-logs/session-YYYY-MM-DD-<topic>.md` describing what changed,
  why, any deviations from spec, open questions, and the state of
  `devtools::test()` and `devtools::check()`. The summary is written for
  review on a train or in bed, not for someone sitting in front of the
  diff. It is staged with the code for Carly's review.
- No session finishes with failing tests or a new check warning. If a
  failure cannot be resolved within the authorized scope, Claude stops,
  documents the problem and its diagnosis in the summary, and surfaces it
  for Carly rather than expanding scope.

## Completed sessions

### Session 1 (completed 58c7e2b): Fix two bugs in null_distribution_F

Infinite loop in the exclude top-up and cryptic base-R error on
`n_random > genome_length`. Both fixed; a single sampling path now handles
truncation with a warning. 15 new tests, zero failures, zero check errors.
NEWS.md entry under "development version, pre-0.1.0".

## Upcoming sessions

### Session 2: Methodological fixes and null-distribution coverage

Authorized changes based on the methodological observations review:

- Item 5: `null_distribution_F` returns the fit-failure count as an
  attribute on the returned vector. The attribute flows through
  `test_decay_significance` and is shown by `print.decay_significance`.
- Item 6: `test_decay_slopes` returns `lme4::isSingular(fit)` and the
  model's convergence code in the `decay_test` object.
  `print.decay_test` prints a one-line diagnostic when the fit is
  singular or non-converged.
- Item 7: `test_decay_slopes` errors (not warns) when any group has fewer
  than two samples, naming the offending group(s) and their counts.
- Item 9: `test_decay_significance` default `exclude_radius` becomes
  `if (is.null(window)) 0 else window`. Behavior change documented in
  roxygen with the recovery option (pass `exclude_radius = 0` explicitly).
- Item 3 (code portion): `null_distribution_F` roxygen revised to
  describe without-replacement sampling accurately. Phipson-Smyth citation
  retained but reframed as the conservative convention, not a derivation.
- Raise `R/null_distribution.R` test coverage above 80%, locking in all
  of the above.

Documentation drafts produced but NOT committed to vignettes:

- Item 3 (vignette portion): draft rewrite of the Phipson-Smyth paragraph
  in `vignettes/null-distribution.Rmd` saved to
  `session-logs/session-YYYY-MM-DD-vignette-drafts.md`.
- Item 4: draft rewrite of the "upper-right cell cannot occur by
  construction" claim in the interpretation table, saved to the same
  file.

Carly edits both drafts substantively and lands them in a separate commit.

Protected-file authorization for Session 2: `R/null_distribution.R`,
`R/decay_tests.R`, and their test files. NOT
`vignettes/null-distribution.Rmd` (drafts only).

Success criteria: all tests pass; check clean; `R/null_distribution.R`
coverage above 80%; each behavior change has a NEWS.md entry; vignette-
draft file populated.

Estimate: one session.

### Session 3: plot_sample fixes, vignette documentation, and palette defaults

Combines Units 2, 3, and 5 from the pre-review plan into one plot-focused
session.

- Fix `plot_sample` defects documented in the earlier plan review:
  x-axis mislabel (plots `bin_mid` but labels "distance from anchor"),
  subtitle inaccuracy, `palette` argument no-op, missing `window`
  default, `log_scale` recomputing instead of using `log_norm_cov`,
  no `sample` validation, flipped argument order
  (`plot_sample(combined, sample, ...)` for API consistency).
- Add "Inspecting a single sample" subsection to `vignettes/intro.Rmd`.
- Switch plot defaults in `plot_coverage_decay`, `plot_decay_slopes`,
  and `plot_null_distribution` to Okabe-Ito (categorical) and Okabe-Ito
  vermillion for the observed-F line in `plot_null_distribution`.
  Preserve existing Brewer palette names for back-compat.

Protected-file authorization: none. This session touches `R/decay_plots.R`,
`tests/testthat/test-decay_plots.R`, `vignettes/intro.Rmd`, `man/`, and
`NEWS.md`.

Success criteria: new tests pass; vignette builds; figures use Okabe-Ito;
check clean.

Estimate: one session.

### Session 4: API surface reduction and minor release hygiene

- Reclassify four helpers as internal: `read_coverage_tab`,
  `bin_coverage`, `normalize_coverage`, `load_sample_coverage`. Remove
  `@export`, add `@keywords internal`, drop `@examples`.
- Keep exported: `add_distance_from_anchor` (user-facing extension point),
  `empirical_pvalue` (needed by the simulation study).
- Regenerate `NAMESPACE` and `man/` via `document()`.
- `.Rbuildignore`: add `simulations`, `NOTES.md`, `PLAN.md`, `CLAUDE.md`,
  `session-logs`, `logs`, `.claude`.
- Remove `viral_rpm` from `globalVariables` (locate it; probably in
  `R/CovDecayAnalyzer-package.R`).

Protected-file authorization: none expected (confirm `globalVariables` is
not in a protected file; route accordingly if it is).

Deferred to a later session: `DESCRIPTION` text corrections (stale ANOVA,
KW, permutation claims; stale RPM reference; `pbkrtest` Imports/Suggests
decision). These require explicit `DESCRIPTION`-edit authorization and
are not in Session 4 scope.

Success criteria: `NAMESPACE` lists only agreed exports; check clean,
no notes about non-standard top-level files; coverage unchanged or higher.

Estimate: one session.

## After these three sessions

Track A2 (simulation study) opens with its own planning document. The
following items depend on simulation output and are scheduled after:

- Kruskal-Wallis plus Dunn benchmark function. Implementation could happen
  earlier but naming and output contract are better fixed with the
  simulation design in view.
- Sample size calculator function, fitted from the simulation power surface.
- `DESCRIPTION` text cleanup (easier once scope is fully settled).

## Open methodological items parked

- Item 8 (residual autocorrelation across adjacent bins; log-transform
  outlier risk from zero-coverage bins with the 1e-5 offset): parked as a
  simulation-design requirement. See `NOTES.md`.
- Item 10 (set.seed state leak in `null_distribution_F`; multiplicity across
  natural anchors; F invariance to distance rescaling): minor polish for a
  future session.
