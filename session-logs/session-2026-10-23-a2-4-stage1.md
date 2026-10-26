# Session summary: A2-4 Stage 1 (Type I production), complete

(Restored: this was first written as `session-2026-10-23-a2-4-interrupted.md`
when Stage 1 finished and Stages 2/3 were held for a budget review; the driver's
stub for the Stage 2 stop later overwrote that file. Content below is the
original Stage 1 summary; the budget-review section is superseded by
`session-2026-10-23-a2-4-engineering.md` and the scope change.)

**Status:** Stage 1 ran to completion with no hard, soft or panic flags. Nothing
committed, no branches, no package code touched.

## Stage 1 results (27 cells, R = 5000, 135,000 datasets, 1,620,000 arm rows)

- 0 failed fits in any arm and cell. Peak worker memory 1.2 GB; cache 0.06 GB.
- **Primary LMM Type I at alpha = 0.05, all 21 null cells: 0.043 to 0.054.**
  Lowest CI lower bound 0.0375, highest CI upper bound 0.0611. Singular rate up
  to 0.54 (lmm), reported per cell in `tbl-typeI-main.csv`.

  | n | G=2 | G=3 | G=4 |
  |---|---|---|---|
  | 2 | 0.0452 | 0.0486 | 0.0522 |
  | 3 | 0.0500 | 0.0544 | 0.0494 |
  | 4 | 0.0534 | 0.0496 | 0.0506 |
  | 5 | 0.0506 | 0.0448 | 0.0444 |
  | 6 | 0.0480 | 0.0478 | 0.0432 |
  | 8 | 0.0542 | 0.0514 | 0.0498 |
  | 10 | 0.0458 | 0.0480 | 0.0430 |

- **Pooled two-stage matches the LMM at every cell** (differences at most 0.005).
- **KW chi-square is off at small n by discreteness:** 0.0 at n = 2, 0.0996 at
  G = 2, n = 3, 0.010 at G = 3, n = 3. KW permutation is 0.0 at G = 2, n <= 3 and
  0.028 to 0.029 at G = 2, n = 4 and 5. The exact chi-square KW null rate is in
  the table beside each KW cell.
- **Post-hoc FWER (`tbl-typeI-posthoc.csv`):** global-null G >= 3: Tukey (LMM)
  0.034 to 0.047 ungated, 0.034 to 0.044 gated; Dunn variants 0.022 to 0.040.
  Partial-null cells (false rejections among the equal pairs): every method at or
  below 0.03.
- Tables: `simulations/results/tier2/typeI/`.

## Compute
Wall 8.82 h (loop version of `kw_perm`; budget then 6.5 h), plus 0.05 h resume.
Six of 1,350 branches errored with "future interrupted" and were re-run by
targets; the driver now retries automatically. The vectorized `kw_perm` returns
identical p-values, so Stage 1 needs no rerun.

## Build notes (A2-4 code, all under `simulations/`)
Sensitivity hooks S1a, S1b, S2, S3; `seed_id` pairing; Tier 2 grids;
`R/sobol.R`; `p_equal_min` column in the arms; `R/summarize_a24.R`;
`R/check_stage_results.R`; `scripts/a24_driver.R`; `tests/testthat/test-a24.R`.
Deviations (collapse cells added because the Sobol design does not vary s_ref;
collapse hard stop requires the CI to exclude 0; S0 has 30 cells; pooled-beats-LMM
flag at effect cells only; S2 as-is analysis only) are in
`session-2026-10-23-a2-4-engineering.md` and the interrupted summary.
