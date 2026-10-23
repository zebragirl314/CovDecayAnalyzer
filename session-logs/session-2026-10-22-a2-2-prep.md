# A2-2 prep (about 5 min of the 15 used)

## 1. Reparameterization: DONE
- `s_i = s_ref * (1 + Delta_std * treated + u_i)`, `u_i ~ N(0, cv^2)`; defaults `Delta_std = 0`, `cv = 1/3` (tau = 2e-6 at s_ref = 6e-6). Edited `R/generate.R`, `R/scenarios.R` (and S3 docstring in `R/sensitivities.R`).
- `make_scenario()` still accepts legacy `Delta=`/`tau=` (per bp) and converts them (`Delta / s_ref`, `tau / s_ref`), so `checks.R` and `pilot_grid()` still run unchanged. Scenario IDs changed (new id fields), so A2-1 pilot IDs no longer match. `pilot_grid()`'s "s_ref x2" collapse cells now mean something different (Delta, tau fixed in per bp at 2x s_ref = half the relative effect); the pilot is not to be rerun.
- `PLAN.md` 2.2, 2.5, 4.1, 4.2 updated. There is no section 4.3; noted in the plan.
- Check: primary LMM, G=2, n=5, R=20, null: Type I 0.05 (CP 0.00-0.25), 0 failures.

## 2. Arms: all four WRITTEN (`analyze_lmm.R` untouched)
`analyze_lmm_slope_only.R`, `analyze_kw_chisq.R`, `analyze_kw_exact.R`, `analyze_two_stage.R`, plus shared `arm_helpers.R`. Same column shape as `analyze_lmm()`.
- KW arms return 3 rows (`arm` = `kw_*`, `kw_*_dunn_ungated`, `kw_*_dunn_gated`). Dunn written by hand (normal approx, Holm); `p_value` on Dunn rows is the minimum Holm-adjusted pairwise p (per-pair p not stored).
- Exact KW enumerates all rank assignments up to 2e5, else 2e4 random permutations (not exact: G=2 n>=9, G=3 n>=5).
- Two-stage: Welch t (G=2, `F_stat` = t^2) or classical one-way ANOVA (G>=3).
- Slope-only LMM fitted directly with lmerTest (no emmeans), same scaling and KR test.

## 3. Dry run (`scripts/dry_run_arms.R`; G=2, n=3, Delta_std=0, 5 kb, R=60, same datasets; `results/explore/a2-2-prep/`)

| Arm | Median s/fit | Errors | Null Type I (CP 95%) | In [0, 0.15] |
|---|---|---|---|---|
| lmm (primary) | 0.139 | 0/60 | 0.10 (0.04-0.21) | yes |
| lmm_slope_only | 0.043 | 0/60 | **0.53 (0.40-0.66)** | **NO** |
| kw_chisq | 0.008 | 0/60 | 0.20 (0.11-0.32) | NO |
| kw_chisq_dunn_ungated / gated | 0.008 | 0/60 | 0.20 / 0.20 | NO |
| kw_exact | 0.008 | 0/60 | 0.00 (0.00-0.06) | yes |
| kw_exact_dunn_ungated / gated | 0.008 | 0/60 | 0.20 / 0.00 | NO / yes |
| two_stage | 0.005 | 0/60 | 0.083 (0.03-0.18) | yes |

KW times include `fit_decay_slopes`; the three KW rows of one fit share one time. No arm errored (flag threshold 5%).

## Recommendation: NO-GO for tonight as things stand
1. **Blocker, wiring.** `run_batch.R` and `_targets.R` still call only `analyze_lmm()`. The new arms are not in the pipeline (outside this session's three tasks). Wiring is small but untested at batch scale.
2. **Slope-only LMM is badly anti-conservative (0.53).** Expected in hindsight: median normalization leaves real between-sample intercept variance (sd_intercept ~ tau * median distance), and dropping the random intercept ignores it. Not a bug I can see; it is a result. It answers the M1 comparator question (the random intercept is needed); drop it from the overnight run or keep as a labelled cautionary arm. Your call.
3. **KW chi-square 0.20 vs 0.10 predicted by M4** at n=3, G=2 (12/60; chance of >=12 under 0.10 is about 2%). Plausibly noise, not verified. `kw_exact` never rejects, as predicted. Worth a look at larger R before trusting that arm.
4. **Dunn-ungated under the exact arm is 0.20** because Dunn uses the normal approximation (G=2 Dunn equals chi-square KW). Only the gate is exact.
5. **1 kb cost** from A2-1 (decision 2) is still unresolved; keep tonight's run at 5 kb.
6. Dry run is one design at R=60; no G>=3 check of the KW/ANOVA arms yet, and `kw_exact` MC fallback is untested.

Not run: anything beyond the dry run; no commits.
