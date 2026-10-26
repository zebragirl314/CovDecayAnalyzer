# Session summary: A2-4 Stage 2 complete, HARD STOP on the collapse check; Stage 3 not started

**Status:** Stage 2 (power) ran to completion (475 cells, 207,000 datasets,
1,164,000 arm rows, 0 errored targets, 0 failed fits). The between-stage check
fired **two collapse hard stops**, so per protocol I halted. Stage 3 was not
started (it was not to auto-start anyway). No panic flags, no soft flags. Nothing
committed, no branches, no package code touched. The Stage 1 summary is in
`session-2026-10-23-a2-4-stage1.md` (it had been overwritten by the driver's
stub; the driver now appends instead).

## The hard stop (needs your decision)

Primary LMM, G = 2, cv = 1/3, Delta_std matched, power across `s_ref`:

| n | Delta_std | s_ref pair | power a | power b | b minus a | 95% CI | flag |
|---|---|---|---|---|---|---|---|
| 5 | 0.5 | 3e-6 vs 6e-6 | 0.5095 | 0.5630 | +0.0535 | +0.016 to +0.091 | HARD (above 0.03, CI excludes 0) |
| 10 | 0.5 | 3e-6 vs 6e-6 | 0.8625 | 0.8970 | +0.0345 | +0.010 to +0.058 | HARD |
| 10 | 0.5 | 3e-6 vs 1.2e-5 | 0.863 | 0.885 | +0.0225 | +0.002 to +0.043 | none (below 0.03) |

The other 15 of 18 comparisons are within 0.03; all Delta_std = 1 comparisons
agree to 0.01 or better (power is near 1). R = 2000 at the two non-reference
`s_ref` levels, R = 1000 at 6e-6. Table: `results/tier2/power/tbl-collapse-tier2.csv`.

### Memo: why I think this is expected, not a fault (half page)

**Concern.** The check as written ("power at matched Delta_std / cv agrees across
s_ref") tests a stronger statement than PLAN 2.5 makes. PLAN 2.5 says power
collapses onto `Delta / tau_eff`, with
`tau_eff = sqrt(tau^2 + sigma_bin^2 / Sxx)`. The second term is within-sample
noise and does **not** scale with `s_ref`, so at fixed Delta_std and cv it
matters more when `s_ref` is small: `tau_eff / tau` is 1.04 at `s_ref` 3e-6,
1.01 at 6e-6 and 1.00 at 1.2e-5.

**Evidence.** A noncentral-F prediction with that `tau_eff` (ncp =
`(Delta/tau_eff)^2 n/2`, df 1 and 2n-2) gives power at n = 5, Delta_std 0.5 of
0.519 (3e-6), 0.542 (6e-6), 0.547 (1.2e-5), and at n = 10 of 0.863, 0.881, 0.885.
Predicted gaps (6e-6 minus 3e-6): +0.023 at n = 5 and +0.018 at n = 10. Observed:
+0.054 and +0.035; these are 1.6 and 1.3 standard errors above the predictions,
so not distinguishable from them. The direction (lower power at the smaller
`s_ref`) matches in every flagged cell.

**Proposed alternative.** Judge the collapse against the `tau_eff` prediction, not
against zero: a hard stop when the observed gap differs from the noncentral-F
predicted gap by more than 0.03 with a CI that excludes the prediction. On that
rule none of the 18 comparisons stops. I did not change the check or the rule; the
stop stands as the prompt defined it, plus my CI qualifier.

**Implications.** (1) The one-input calculator should be indexed by
`Delta / tau_eff`, not by Delta_std alone; ignoring the noise term costs up to
about 0.05 in power at the smallest `s_ref` (3e-6) and under 0.02 at the larger
two, at 5 kb bins with the reference sigma_atom 0.15. (2) NOTES.md's "Collapse
hypothesis supported at Tier 1... indexed by the standardized effect alone" is
slightly too strong; Tier 1 at R = 200 could not see a 0.02 effect. This is for
you to decide; I changed no methodological text. (3) The surface fit (a later
session) tests this directly with the `Delta / tau_eff` covariate.

## Stage 2 headline numbers

- **0 failed fits in any arm and cell.** Peak worker memory 1.43 GB, cache 0.06
  GB, wall 7.51 h of the 10 h budget (15 h panic threshold); three branches were
  "future interrupted" again and were retried automatically.
- **Primary LMM power, G = 2, Delta_std = 1 (reference cv):** n = 2: 0.381,
  n = 3: 0.769, n = 5: **0.985**, n = 10: 1.000. G = 3: 0.441, 0.821, 0.987;
  G = 4: 0.467, 0.841, 0.995 (n = 2, 3, 5).
- **Pooled two-stage vs LMM:** identical decisions within 0.006 on all 63 figure
  cells (paired difference pooled minus LMM from -0.002 to +0.006, mean about 0).
  No "pooled beats LMM" soft flag. Welch two-stage at G = 2 loses at small n
  (n = 2: 0.186 vs 0.382; n = 3: 0.665 vs 0.769), as in A2-3.
- **KW:** chi-square and permutation have zero or erratic power at n = 2 and 3
  (for example G = 2, n = 3, Delta_std 1: chi-square 0.880, permutation 0.000;
  G = 3, n = 3: 0.164 and 0.579), by discreteness.
- **Smallest n for 80% power, LMM, reference cv (`tbl-n-for-80pct.csv`):**

  | G | Delta_std 0.5 | 1.0 | 1.5 |
  |---|---|---|---|
  | 2 | 10 (power 0.897) | 4 | 3 |
  | 3 | 8 | 3 (CI-conservative: 4) | 3 |
  | 4 | 8 | 3 | 3 |

- **Sobol surface (`surface/power-surface-data.csv`):** 400 points, LMM only, 0
  failures, mean singular rate 0.16.
- Files in `simulations/results/tier2/power/`: `tbl-power-main.csv`,
  `tbl-n-for-80pct.csv`, `tbl-collapse-tier2.csv`, `grid.csv`,
  `surface/power-surface-data.csv`, `raw/tier2-power.rds` (xz).

## Stage 3 status and projection (not run)

Configured as instructed: R = 1000 for all sensitivity cells (222 cells, 222,000
datasets), slope-only arm not called (11 arm/variant rows per dataset; documented
in the `a24_sens_grid()` docstring), budget 6 h (panic 9 h), session cap 24 h.
Projection from the earlier throughput sample (0.68 core-s per dataset with six
arms, about 0.56 without the slope-only arm): about 3.5 h, range 2.5 to 5 h; fits.
Resume with `caffeinate -d -i -s Rscript scripts/a24_driver.R sensitivities`.
Note the driver stops on any Stage 3 hard flag (Type I CI above 0.08 or below
0.025, an arm failing in more than 10% of fits), which genuine sensitivity
effects such as dropouts or heteroscedasticity can trigger; if that happens the
tables are already written before the check runs.

## Open questions

- Accept the `tau_eff`-based reading of the collapse check (memo above) and treat
  the two hard stops as resolved? Or keep the stop as a finding for the paper?
- Edit the NOTES.md collapse wording? (Protected section; not touched.)
- Fix for the process: the driver wrote its stub over an existing summary; it now
  appends. Dropping the file name collision entirely would mean per-stage names;
  say if you want that.
# A2-4 interrupted (stub written by the driver)

Stage: sensitivities

Reason: HARD STOP

### Stage status: sensitivities

- Cells run: 222
- Wall time: 2.67 h
- Continue: FALSE
- Hard-stop flags fired (5):
  - lmm_typeI_ci_below_0.025: sens_S2 G=2 n=3 Delta_std=0 [dropout_frac_runkb=0.02_5] (value 0.01831; CI upper < 0.025). rate 0.0100, CP 0.0048 to 0.0183, R = 1000
  - lmm_typeI_ci_below_0.025: sens_S2 G=3 n=3 Delta_std=0 [dropout_frac_runkb=0.02_5] (value 0.01163; CI upper < 0.025). rate 0.0050, CP 0.0016 to 0.0116, R = 1000
  - lmm_typeI_ci_below_0.025: sens_S2 G=2 n=5 Delta_std=0 [dropout_frac_runkb=0.02_5] (value 0.0196; CI upper < 0.025). rate 0.0110, CP 0.0055 to 0.0196, R = 1000
  - lmm_typeI_ci_below_0.025: sens_S2 G=2 n=3 Delta_std=0 [dropout_frac_runkb=0.05_5] (value 0.01301; CI upper < 0.025). rate 0.0060, CP 0.0022 to 0.0130, R = 1000
  - lmm_typeI_ci_below_0.025: sens_S2 G=3 n=3 Delta_std=0 [dropout_frac_runkb=0.05_5] (value 0.01437; CI upper < 0.025). rate 0.0070, CP 0.0028 to 0.0144, R = 1000
- Soft-stop flags fired (25):
  - lmm_typeI_ci_straddles_0.07: sens_S4 G=3 n=3 Delta_std=0 [cv_ratio=2] (value 0.067; CI contains 0.07, estimate in [0.06, 0.07]). rate 0.0670, CP 0.0523 to 0.0843, R = 1000
  - lmm_typeI_ci_straddles_0.07: sens_S4 G=3 n=5 Delta_std=0 [cv_ratio=2] (value 0.07; CI contains 0.07, estimate in [0.06, 0.07]). rate 0.0700, CP 0.0550 to 0.0876, R = 1000
  - lmm_typeI_ci_straddles_0.07: sens_S5 G=3 n=3 Delta_std=0 [anchor_frac=1] (value 0.063; CI contains 0.07, estimate in [0.06, 0.07]). rate 0.0630, CP 0.0487 to 0.0799, R = 1000
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=5 Delta_std=0.5 [dropout_frac_runkb=0.005_5] (value 0.034; pooled - LMM > 0.02). pooled 0.3690 vs LMM 0.3350; 34 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=3 Delta_std=1 [dropout_frac_runkb=0.005_5] (value 0.047; pooled - LMM > 0.02). pooled 0.6090 vs LMM 0.5620; 49 pooled-only, 2 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=3 Delta_std=0.5 [dropout_frac_runkb=0.02_5] (value 0.065; pooled - LMM > 0.02). pooled 0.1070 vs LMM 0.0420; 65 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=3 Delta_std=0.5 [dropout_frac_runkb=0.005_5] (value 0.039; pooled - LMM > 0.02). pooled 0.1970 vs LMM 0.1580; 39 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=5 Delta_std=1 [dropout_frac_runkb=0.05_5] (value 0.031; pooled - LMM > 0.02). pooled 0.1480 vs LMM 0.1170; 32 pooled-only, 1 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=3 Delta_std=1 [dropout_frac_runkb=0.005_5] (value 0.059; pooled - LMM > 0.02). pooled 0.6140 vs LMM 0.5550; 59 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=5 Delta_std=1 [dropout_frac_runkb=0.02_5] (value 0.072; pooled - LMM > 0.02). pooled 0.4300 vs LMM 0.3580; 74 pooled-only, 2 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=5 Delta_std=0.5 [dropout_frac_runkb=0.05_5] (value 0.028; pooled - LMM > 0.02). pooled 0.0740 vs LMM 0.0460; 28 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=5 Delta_std=0.5 [dropout_frac_runkb=0.005_5] (value 0.046; pooled - LMM > 0.02). pooled 0.3980 vs LMM 0.3520; 46 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=5 Delta_std=1 [dropout_frac_runkb=0.05_5] (value 0.06; pooled - LMM > 0.02). pooled 0.2140 vs LMM 0.1540; 61 pooled-only, 1 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=3 Delta_std=1 [dropout_frac_runkb=0.05_5] (value 0.08; pooled - LMM > 0.02). pooled 0.1270 vs LMM 0.0470; 80 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=3 Delta_std=1 [dropout_frac_runkb=0.02_5] (value 0.112; pooled - LMM > 0.02). pooled 0.3330 vs LMM 0.2210; 113 pooled-only, 1 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=3 Delta_std=0.5 [dropout_frac_runkb=0.005_5] (value 0.023; pooled - LMM > 0.02). pooled 0.2120 vs LMM 0.1890; 23 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=3 Delta_std=1 [dropout_frac_runkb=0.02_5] (value 0.106; pooled - LMM > 0.02). pooled 0.3030 vs LMM 0.1970; 107 pooled-only, 1 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=3 Delta_std=0.5 [dropout_frac_runkb=0.05_5] (value 0.033; pooled - LMM > 0.02). pooled 0.0460 vs LMM 0.0130; 34 pooled-only, 1 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=5 Delta_std=0.5 [dropout_frac_runkb=0.05_5] (value 0.036; pooled - LMM > 0.02). pooled 0.0760 vs LMM 0.0400; 36 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=3 Delta_std=0.5 [dropout_frac_runkb=0.05_5] (value 0.056; pooled - LMM > 0.02). pooled 0.0730 vs LMM 0.0170; 56 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=3 Delta_std=1 [dropout_frac_runkb=0.05_5] (value 0.078; pooled - LMM > 0.02). pooled 0.1300 vs LMM 0.0520; 78 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=3 n=5 Delta_std=0.5 [dropout_frac_runkb=0.02_5] (value 0.047; pooled - LMM > 0.02). pooled 0.1250 vs LMM 0.0780; 49 pooled-only, 2 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=3 Delta_std=0.5 [dropout_frac_runkb=0.02_5] (value 0.046; pooled - LMM > 0.02). pooled 0.1170 vs LMM 0.0710; 46 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=5 Delta_std=0.5 [dropout_frac_runkb=0.02_5] (value 0.059; pooled - LMM > 0.02). pooled 0.1670 vs LMM 0.1080; 60 pooled-only, 1 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S2 G=2 n=5 Delta_std=1 [dropout_frac_runkb=0.02_5] (value 0.081; pooled - LMM > 0.02). pooled 0.4830 vs LMM 0.4020; 82 pooled-only, 1 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
- Panic-exit flags fired (0):
  none


See `session-logs/session-2026-10-23-a2-4.md` for the status blocks of the stages that completed, and `a2-4-resource-trace.csv`.

---
## Stage 3 stop (2026-10-24, 5 S2 hard stops): reclassified 2026-10-25
See the Stage 3 section of `session-2026-10-23-a2-4.md`. All three stages are
complete. Nothing committed.
