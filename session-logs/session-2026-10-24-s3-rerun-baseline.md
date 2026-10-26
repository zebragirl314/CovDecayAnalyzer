# Session summary: S3 controlled baseline (sigma_ratio = 1) at sigma_atom 0.3 and 0.5

**Status:** complete, clean. 24 cells (G in {2, 3} x n in {3, 5} x Delta_std in
{0, 0.5, 1} x sigma_atom in {0.3, 0.5}), homoscedastic, R = 1000, 24,000 datasets,
264,000 rows, 0 errored targets, 0 failed fits. Wall 0.25 h (budget 0.5 h). Peak
worker memory 0.91 GB, cache 0.007 GB. **0 hard, 0 soft, 0 panic flags.** Nothing
committed. Edits: grid `a24_s3baseline_grid()` (`R/scenarios.R`), writer
`write_a24_s3baseline()` (`R/summarize_a24.R`), stage `s3baseline`
(`scripts/a24_driver.R`), and the annotation on `R/check_stage_results.R` marking the
cells as the controlled baseline. Output:
`results/tier2/sensitivities/tbl-sens-hetero-residual-baseline.csv` (LMM, pooled,
Welch rate with CP CI; paired LMM minus pooled; LMM minus Welch),
`grid-s3baseline.csv`, `raw/tier2-s3baseline.rds`. The baseline cells share draws
with the sigma_ratio = 2 rerun cells (scenario_id = their seed_id).

## Pre-committed interpretation, applied

Rule: agree at sigma_ratio 1 while differing at 2 means heteroscedasticity; still
differing at 1 means noise level. **Result: mostly the first, with one residual
case of the second.** LMM minus pooled has a 95% CI excluding 0 in **2 of 24**
cells at sigma_ratio 1, against **10 of 24** at sigma_ratio 2. Largest gap -0.015
at ratio 1, -0.041 at ratio 2. The two baseline cells are both G=2, n=3,
sigma_atom 0.5 (Delta_std 0: -0.007, CI -0.014 to -0.001; Delta_std 0.5: -0.015,
CI -0.023 to -0.007). I am not rounding this to either branch: the divergence is
mainly a heteroscedasticity effect, and there is a smaller noise-level effect at the
smallest design at the largest noise.

## Headline sentences

1. **At matched sigma_atom with homoscedastic residuals, LMM and pooled agree
   almost everywhere (within 0.005 in 22 of 24 cells; G=3 cells all within 0.004),
   and diverge only at G=2, n=3, sigma_atom 0.5**, where the pooled test rejects
   more by 0.007 (Delta_std 0) and 0.015 (Delta_std 0.5), against 0.025 and
   0.041 with heteroscedastic residuals in the same cells. At G=3, n=3, sigma 0.5
   the gap is -0.004 and -0.003 (CIs cover 0) at ratio 1 against -0.015 and
   -0.024 (CIs exclude 0) at ratio 2.
2. **The LMM is mildly below nominal at ratio 1, not severely conservative.**
   Type I is 0.029 to 0.045 (mean 0.037); pooled 0.034 to 0.049. At ratio 2
   the LMM range was 0.017 to 0.046. The lowest ratio-1 cell is 0.029 (CP 0.020
   to 0.041), G=2, n=3, sigma 0.5; four cells (that one, and G=2 n=5, G=3 n=3 and G=3 n=5 at
   sigma 0.3) have CIs entirely below 0.05. No cell is near the 0.025
   stop.
3. **Large noise alone makes the LMM fit singular but, in these designs, not
   different from the pooled test.** At Delta_std 0 the singular-fit rate is
   0.19 to 0.50 at sigma 0.3 and 0.66 to 0.78 at sigma 0.5 (Stage 3 reference
   at 0.15: 0.01 to 0.21), yet LMM and pooled rates differ by at most 0.007 at
   Type I cells. Residual-SD ratio 2 on top raises the median fitted residual SD
   (0.35 to 0.40 vs 0.23 at sigma 0.5) and the singular rate (0.87 to 0.92) and is
   where the LMM falls behind.
4. **Power cost of each factor, LMM, G=2, n=3, Delta_std 1:** noise alone: 0.769
   (sigma 0.15, Stage 3, unpaired) to 0.731 (0.3) to 0.695 (0.5); adding ratio 2:
   0.688 and 0.569. Pooled loses about the same amount to ratio 2 (0.730 to
   0.689, 0.697 to 0.586), so the ratio-2 power loss itself is shared by both tests;
   the LMM-specific part is the extra gap of 0.002 (sigma 0.3) to 0.015 (sigma 0.5) in
   these cells (-0.001 vs +0.001; -0.017 vs -0.002).
5. **Welch two-stage at G=2** trails the LMM at n=3 by 0.056 and 0.104 (Delta_std
   0.5, 1, sigma 0.3) and 0.049 and 0.113 (sigma 0.5), as in the balanced reference
   (degrees of freedom at tiny n); at G=3 it is the pooled test.

## What this means for the paper

- Heteroscedastic residuals (SD ratio 2) are what separates the LMM from the pooled
  test at n=3, and the separation grows with the noise level (none at
  sigma_atom 0.15; up to 0.015 at 0.3; 0.014 to 0.041 at 0.5, all n=3 cells). It runs in the pooled
  test's favour. The recommendation is about heteroscedastic, noisy, very
  small designs (n=3), not about noise level in general.
- At n=5 the divergence is at most 0.013 with ratio 2 (G=2, sigma 0.5, Delta_std
  0.5) and at most 0.003 at ratio 1.
- One rerun, one interpretation; no further S3 variations without reopening the
  question.

## Caveats
- The ratio-1 versus ratio-2 contrasts per cell are on the same draws; the
  versus-Stage-3-reference (sigma_atom 0.15) contrasts are unpaired.
- The residual ratio-1 gaps (-0.015, CI -0.023 to -0.007; -0.007, CI -0.014 to -0.001)
  are small and sit in one design (G=2, n=3, sigma 0.5), both in the same direction;
  two of 24 intervals excluding 0 is about what chance alone gives, so treat
  them as suggestive.
- The two S3 reruns changed the last group's residual SD only; the mechanism
  (fits going singular as noise rises) is seen in the fits, not tested.

---

### Stage status: sensitivities

- Cells run: 24
- Wall time: 0.25 h
- Continue: TRUE
- Hard-stop flags fired (0):
  none
- Soft-stop flags fired (0):
  none
- Panic-exit flags fired (0):
  none

Peak worker memory (MB): 912.1; cache size (GB): 0.007; datasets: 24000; files: tbl-sens-hetero-residual-baseline.csv, grid-s3baseline.csv, tier2-s3baseline.rds

