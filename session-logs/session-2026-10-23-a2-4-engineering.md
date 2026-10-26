# Session summary: A2-4 engineering fix (vectorized kw_perm, LMM-only surface points)

**Status:** Both tasks and the targets change are done and tested. Stages 2 and 3
were NOT resumed (waiting for your review of the timing below). Nothing committed,
no branches, no package code touched. Writes confined to `simulations/` and
`session-logs/`. `simulations/tests/testthat` passes in full (a24, arms, generator,
kw-perm-invariant). `devtools::test()`/`check()` not run (no package changes).

## Task 1: vectorized random-permutation path (`R/analyze_kw_perm.R`)

- `kw_perm_p()` random-permutation branch now builds one N x B label matrix and
  takes column sums per group, instead of `replicate(B, tapply(...))`. The omnibus
  statistic over the permutations is also vectorized (`rowSums`). The permutation
  Dunn statistics already derive from the same `sums` matrix, so omnibus and
  Dunn-perm both use the vectorized path.
- The permutations are drawn with the same calls in the same order as before
  (`sample(g)` is `g[sample.int(N)]`), so the RNG stream is unchanged and the new
  path returns **identical** p-values, not only the same distribution. The
  exact-enumeration branch is untouched (it keeps its `apply` statistic). I used
  base R; `matrixStats` is not in the renv library and adding it is a dependency
  decision.
- One-dataset check, G = 3, n = 5: 0.058 s vs 0.507 s (about 9x); `identical()` TRUE.
- **Consequence: Stage 1 results are not affected.** Stage 1 ran the loop version;
  the vectorized version reproduces its p-values exactly, so no Stage 1 rerun is
  needed and A2-2/A2-3 results stay reproducible.

**Dormant reference for reproducibility** (same file): `kw_perm_p_loop_internal()`
(the loop implementation, verbatim) and `analyze_kw_perm_loop_internal()` (the
wrapper). Not used by the pipeline.

## Task 2: invariant test (`tests/testthat/test-kw-perm-invariant.R`)

G = 3, n = 5, Delta_std = 0, reference nuisance, 500 datasets (replicates 1 to
500), every fit on the random-permutation path (verified), none failed.

| Variant | vectorized | loop | CP overlap | max abs p diff |
|---|---|---|---|---|
| omnibus | 26/500 (CP 0.034 to 0.075) | 26/500 (same) | 1.00 | 0 |
| dunn_perm_ungated | 23/500 (0.029 to 0.068) | 23/500 (same) | 1.00 | 0 |
| dunn_perm_gated | 23/500 (0.029 to 0.068) | 23/500 (same) | 1.00 | 0 |

**Passes.** The requested criterion (CP intervals overlap by at least 50%) is
asserted, plus a stronger one (p-values equal), which holds because the streams are
the same. A second test compares the two functions directly at G = 4, n = 4 and
confirms the enumeration path (G = 2, n = 3) is identical. Note that because the
streams coincide the 50%-overlap criterion is met trivially; it cannot detect a
change in the permutation distribution, only the equality assertion would.

## Targets change: Sobol surface points run the LMM only

- `run_batch()` gained an `arms` argument (default: all six, so every existing call
  is unchanged); `run_batch_slim()` passes it through; `a24_arms_for(batch)` returns
  `"lmm"` for `power_surface` batches and all arms otherwise; `_targets.R` calls it.
  The 63 figure cells, the 12 collapse cells and all sensitivity families keep all
  arms (12 arm/variant rows per dataset).
- Follow-on edits so the driver and checks still work: expected row count in
  `scripts/a24_driver.R` is 12 rows per dataset, 1 for surface points; the
  "pooled beats LMM" check only looks at cells where both arms ran. Test added in
  `test-a24.R`. Note: "seven variants" in your message is 12 arm/variant rows
  (six arm families); surface points now produce 1 row per dataset.
- One-cell dry run (G = 3, n = 3, R = 10, reference): 120 rows, 0 failed, all
  12 arm/variant rows present with 10 finite p-values each.

## New timing estimates

Measured with the vectorized path (A2-4 `run_batch_slim`, correct arms per cell,
BLAS pinned to 1 thread), machine otherwise idle. You asked for three to five
sample fits; I used more because a handful of cells is dominated by which cells
happen to be drawn (cells sampled in proportion to replicate count).

| Measurement | Stage 1 (loop, check) | Stage 1 (vectorized) | Stage 2 | Stage 3 |
|---|---|---|---|---|
| Serial, mean s per dataset (27 / 27 / 40 / 40 cells x 2 reps) | 0.735 | 0.425 | 0.403 | 0.287 |
| 10 processes, core-s per dataset (200 to 300 datasets) | 1.42 | 1.11 | 1.26 | 0.68 |
| Projected wall at 10 workers | 5.3 h | 4.2 h | 7.3 h | 8.4 h |
| Sum of per-arm times, projected wall | | | 4.9 h | 5.5 h |
| Calibrated (x1.66, see below) | 8.8 h (actual) | 7.0 h | about 12 h | about 14 h |

Calibration: the real Stage 1 (loop) took 8.82 h against 5.3 h projected from the
same kind of sample, so short samples understate by about 1.66x (targets overhead,
stragglers, other load on the Mac). Read "projected" as the optimistic end and
"calibrated" as the pessimistic end. Central estimates: **Stage 2 about 7 h
(range 5 to 12 h); Stage 3 about 8 h (range 5.5 to 14 h).**

Where the time is now (10-process sample, s per dataset): Stage 2 LMM 0.75
(surface points 0.95, figure cells 0.53, collapse cells 0.41), `kw_perm` 0.12,
slope-only LMM 0.06, everything else under 0.02. Stage 3 LMM 0.29, slope-only LMM
0.08, `kw_perm` 0.04, rest under 0.02. **The primary LMM fit is now 70 to 80% of
the cost**; the vectorization removed the `kw_perm` hotspot but the remainder
cannot be sped up without touching `R/decay_tests.R`.

## Do the panic budgets fit?

**No.** Panic thresholds are 150% of budget: Stage 2 2.25 h, Stage 3 4.8 h, 20 h
for the session. Even the optimistic estimates (Stage 2 4.9 h, Stage 3 5.5 h)
exceed them, and Stage 1 already used 8.9 h, so the 20 h session total would also
be crossed (8.9 + 7 + 8 = about 24 h central). Running as specified would trigger
a panic exit in both stages.

What would fit (hours scale linearly in datasets; central estimates):

- Stage 3 at R = 1000 instead of 2000: about 4 h central (2.8 to 7 h), still near
  the 4.8 h threshold. Dropping the slope-only comparator from Stage 3 saves about
  18% of the cost.
- Stage 2: surface points are 120,000 of the 207,000 datasets and about 3.2 h at
  0.95 s each. Surface R = 150 and collapse cells R = 1000 would make Stage 2 about
  3.5 to 4 h.
- Or raise the budgets: a realistic ask is Stage 2 about 10 h cap, Stage 3 about 14 h
  cap, session cap 30 h.

These are decisions for you; I changed none of the budgets or replicate counts.

## Files changed

`R/analyze_kw_perm.R`, `R/run_batch.R`, `R/summarize_a24.R`, `R/check_stage_results.R`,
`_targets.R`, `scripts/a24_driver.R`, `tests/testthat/test-a24.R`,
new `tests/testthat/test-kw-perm-invariant.R`. Timing scripts are in the session
scratchpad, not the repo.

## Open questions

- Which budget or replicate-count option for Stages 2 and 3?
- Keep the slope-only cautionary arm in Stage 3 (about 18% of its cost)?
- Resume command once decided (from `simulations/`):
  `caffeinate -dims Rscript scripts/a24_driver.R power sensitivities`
  (edit `STAGES` budgets and the `R` arguments in `a24_power_grid()` and
  `a24_sens_grid()` first if you change them).
