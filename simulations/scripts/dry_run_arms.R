## A2-2 prep dry run. Run from simulations/:  Rscript scripts/dry_run_arms.R
## (1) Reparameterization check: primary LMM, reference geometry, G=2, n=5,
##     null, R=20 (the pilot's cell).
## (2) All five arms on R=60 replicates of the smallest grid scenario
##     (n=3, G=2, Delta_std=0, 5 kb bins, reference geometry), same datasets.
source("R/setup.R")
load_sim()
out_dir <- "results/explore/a2-2-prep"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

typeI <- function(p, n) {
    rej <- sum(p <= 0.05, na.rm = TRUE)
    ci <- cp_interval(rej, n)
    c(rate = rej / n, lower = ci[["lower"]], upper = ci[["upper"]])
}

sc5 <- make_scenario(G = 2L, n = 5L, Delta_std = 0)
chk <- do.call(rbind, lapply(1:20, function(r) analyze_lmm(generate_replicate(sc5, r))))
t5 <- typeI(chk$p_value, 20)
cat(sprintf("Reparam check (primary, G=2 n=5, R=20): Type I %.2f (CP %.2f-%.2f), failed %d\n",
    t5[["rate"]], t5[["lower"]], t5[["upper"]], sum(chk$status != "ok")))

sc <- make_scenario(G = 2L, n = 3L, Delta_std = 0)
arms <- list(lmm = analyze_lmm, lmm_slope_only = analyze_lmm_slope_only,
    kw_chisq = analyze_kw_chisq, kw_exact = analyze_kw_exact,
    two_stage = analyze_two_stage)
invisible(lapply(arms, function(f) f(generate_replicate(sc, 999L))))  # warm-up
R <- 60L
res <- do.call(rbind, lapply(seq_len(R), function(r) {
    d <- generate_replicate(sc, r)
    do.call(rbind, lapply(arms, function(f) f(d)))
}))
utils::write.csv(res, file.path(out_dir, "dry-run-replicates.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(res, res$arm), function(d) {
    ti <- typeI(d$p_value, nrow(d))
    tibble::tibble(arm = d$arm[1L], fits = nrow(d),
        error_rate = mean(d$status != "ok"),
        wall_median_s = stats::median(d$wall_time, na.rm = TRUE),
        wall_p95_s = unname(stats::quantile(d$wall_time, .95, na.rm = TRUE)),
        typeI = ti[["rate"]], ci_lower = ti[["lower"]], ci_upper = ti[["upper"]],
        in_0_0.15 = ti[["rate"]] <= 0.15, flag_errors = mean(d$status != "ok") > 0.05)
}))
print(as.data.frame(summ), digits = 3)
utils::write.csv(summ, file.path(out_dir, "dry-run-summary.csv"), row.names = FALSE)
print(table(res$status[res$status != "ok"]))
