# Session summary: S3 rerun at larger sigma_atom (A2-4 follow-up)

**Status:** complete, clean. 24 cells (G in {2, 3} x n in {3, 5} x Delta_std in
{0, 0.5, 1}, residual SD ratio 2, sigma_atom in {0.3, 0.5}), R = 1000, 24,000
datasets, 264,000 rows (11 per dataset, no slope-only arm), 0 errored targets, 0
failed fits. Wall 0.25 h (budget 1 h; faster than the one hour expected). Peak worker memory 0.86 GB, cache
0.007 GB. **0 hard stops, 0 panic, 2 soft flags** (pooled beats LMM by more than
0.02; below). Nothing committed. The only edit outside the run code is the
annotation on `R/check_stage_results.R` marking these cells as a documented
follow-up on the original S3. Run code: `a24_s3rerun_grid()` in `R/scenarios.R`,
`write_a24_s3rerun()` in `R/summarize_a24.R`, stage `s3rerun` in
`scripts/a24_driver.R` (checked with the "sensitivities" thresholds).
Output: `results/tier2/sensitivities/tbl-sens-hetero-residual-rerun.csv`
(per cell: LMM, pooled, Welch omnibus rate with CP CI; paired LMM minus pooled with
Newcombe CI and discordant counts; LMM minus Welch), `grid-s3rerun.csv`,
`raw/tier2-s3rerun.rds`.

## Headline sentences

1. **Does the noise term reach the test? Yes at sigma_atom 0.5 and at n = 3,
   weakly at 0.3.** tau_eff/tau is only 1.04 (0.3) and 1.10 (0.5), so the
   standardized effect barely changes, but the LMM's behaviour does: median fitted
   residual SD is 0.20 to 0.22 (0.3) and 0.35 to 0.40 (0.5), against 0.07 at
   sigma_atom 0.15, and the singular-fit rate at Delta_std 0 rises from 0.01 to
   0.21 (Stage 3 reference) to 0.53 to 0.78 (0.3) to 0.87 to 0.92 (0.5). LMM power at
   G=2, n=3, Delta_std 1 is 0.688 (0.3) and 0.569 (0.5), against 0.769 in the
   Stage 3 reference (different noise, unpaired, no same-sigma baseline cell was
   run, so the change cannot be attributed to the SD ratio alone).
2. **Do LMM and pooled still agree? At sigma_atom 0.3, yes except at G=2, n=3;
   at 0.5 no, at n = 3: the pooled test rejects more.** At 0.3, LMM minus pooled
   is within 0.003 in 10 of 12 cells; exceptions G=2 n=3, Delta_std 0 (-0.015) and
   0.5 (-0.007), CIs exclude 0. At 0.5, LMM minus pooled is -0.014 to -0.041 in all
   n=3 cells (CIs exclude 0) and -0.009 and -0.013 at G=2 n=5 (Delta_std 0, 0.5);
   it is within 0.003 at G=3 n=5 and at G=2 n=5, Delta_std 1. The largest gaps are
   G=2, n=3, Delta_std 0.5 (0.145 vs 0.186, -0.041, 43 pooled-only vs 2 LMM-only
   rejections) and Delta_std 0 (0.017 vs 0.042, -0.025).
3. **The LMM is conservative, not inflated, here too.** LMM Type I is 0.034 to
   0.046 at 0.3 and 0.017 to 0.046 at 0.5, with the lowest cell (G=2, n=3, sigma
   0.5) at 0.017 (CP 0.010 to 0.027; the upper bound is above 0.025, so no
   hard stop). The pooled test is 0.034 to 0.061 (0.061 at G=3 n=3, sigma 0.5; the
   upper CI is not near a stop). Same pattern as S2 dropouts: the LMM loses power when the
   within-sample noise is large relative to the between-sample slope variance and
   the fit goes singular. Mechanism seen in the fits only; I did not test why the
   pooled test is less affected.
4. **Welch two-stage at G=2:** LMM minus Welch is +0.042 and +0.112 at sigma 0.3, G=2,
   n=3, Delta_std 0.5 and 1 (pooled variance buys degrees of freedom at n=3), +0.080
   at 0.5, n=3, Delta_std 1; at G=3 Welch and pooled are the same test (differences
   equal).

## Soft flags (2)
`pooled_two_stage_beats_lmm_gt_0.02`: G=2 n=3 Delta_std 0.5, sigma 0.5 (pooled 0.186
vs LMM 0.145); G=3 n=3 Delta_std 0.5, sigma 0.5 (0.224 vs 0.200). Both are the
n=3, large-noise condition above; they extend the S2 finding, they do not contradict it.

## What this means for the paper
Heteroscedastic residuals alone (S3 at reference sigma) do not separate the LMM from
the pooled test (Stage 3 section, `session-2026-10-23-a2-4.md`). With enough
within-sample noise (sigma_atom 0.5, n=3) they do, in the pooled test's favour by
0.01 to 0.04 in rejection rate. This is the same direction as S2. No cell in
Stage 3 or this rerun shows the LMM beating the pooled test.

## Caveats
- No sigma_ratio = 1 baseline at sigma 0.3 and 0.5, as scoped; the effect of the
  residual SD ratio itself is not separated from the effect of larger sigma_atom.
  A 24-cell baseline (about 15 min at this rate) would do it; say if you want it.
- Paired differences are between arms on the same datasets; comparisons with
  Stage 3 reference cells are unpaired and across different sigma.
- A fixed 5 kb bin and the reference L, so tau_eff/tau stays near 1; larger
  tau_eff effects need smaller s_ref (Stage 2 collapse cells), not tested here.

---

### Stage status: sensitivities

- Cells run: 24
- Wall time: 0.25 h
- Continue: TRUE
- Hard-stop flags fired (0):
  none
- Soft-stop flags fired (2):
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S3rerun G=2 n=3 Delta_std=0.5 [sigma_atom=0.5] (value 0.041; pooled - LMM > 0.02). pooled 0.1860 vs LMM 0.1450; 43 pooled-only, 2 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
  - pooled_two_stage_beats_lmm_gt_0.02: sens_S3rerun G=3 n=3 Delta_std=0.5 [sigma_atom=0.5] (value 0.024; pooled - LMM > 0.02). pooled 0.2240 vs LMM 0.2000; 24 pooled-only, 0 LMM-only of 1000 (CONTRADICTS the Tier 1 pivot)
- Panic-exit flags fired (0):
  none

Peak worker memory (MB): 856.7; cache size (GB): 0.007; datasets: 24000; files: tbl-sens-hetero-residual-rerun.csv, grid-s3rerun.csv, tier2-s3rerun.rds

