## Tier 0 pilot (PLAN.md section 3). Run from simulations/:
##   Rscript scripts/pilot_timing.R
##
## Serial by design (clean timing). R = 20 per cell unless overridden:
##   CDA_PILOT_R           replicates per cell (default 20)
##   CDA_PILOT_R_M3        replicates for the M3 anchor cells (default: R)
##   CDA_PILOT_R_COLLAPSE  replicates for the collapse cells (default: R)
##   CDA_PILOT_OUT         output directory (default results/explore/pilot)
##
## Measures: (1) wall time per test_decay_slopes call, (2) singular rate (M1),
## (3) failure rate, (4) Type I at three anchor positions (M3), (5) the
## collapse hypothesis at one point (2.5). Also records df2 by design (M2),
## the emmeans share of call time (M6), convergence codes (M7), and the
## arithmetic behind M4. Writes CSVs and pilot-report.md.

source("R/setup.R")
load_sim()

R_default  <- as.integer(Sys.getenv("CDA_PILOT_R", "20"))
R_m3       <- as.integer(Sys.getenv("CDA_PILOT_R_M3", R_default))
R_collapse <- as.integer(Sys.getenv("CDA_PILOT_R_COLLAPSE", R_default))
out_dir    <- Sys.getenv("CDA_PILOT_OUT", "results/explore/pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- pilot_grid()
cell_R <- function(purposes) {
    tags <- strsplit(purposes, ";", fixed = TRUE)[[1L]]
    max(c(
        if (any(tags == "m3")) R_m3,
        if (any(grepl("^collapse", tags))) R_collapse,
        if (any(!tags %in% "m3" & !grepl("^collapse", tags))) R_default
    ))
}
grid$R <- vapply(grid$purposes, cell_R, integer(1L))

## ---- run every cell, serially -----------------------------------------
invisible(analyze_lmm(generate_replicate(make_scenario(), 1L)))  # warm-up

t_start <- proc.time()[["elapsed"]]
res_list <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    batch <- cbind(row, rep_start = 1L, rep_end = row$R)
    t0 <- proc.time()[["elapsed"]]
    res_list[[i]] <- run_batch(batch)
    res_list[[i]]$purposes <- row$purposes
    cat(sprintf("[%2d/%d] %-30s G=%d n=%2d L=%6.0f bin=%4d R=%3d  %6.1f s\n",
        i, nrow(grid), row$purposes, row$G, row$n, row$L, row$bin_size,
        row$R, proc.time()[["elapsed"]] - t0))
}
sim <- do.call(rbind, res_list)
cat(sprintf("Total pilot time: %.1f s over %d fits\n",
    proc.time()[["elapsed"]] - t_start, nrow(sim)))
utils::write.csv(sim, file.path(out_dir, "pilot-replicates.csv"),
    row.names = FALSE)

has_tag <- function(d, tag) {
    vapply(strsplit(d$purposes, ";", fixed = TRUE), function(t) any(t == tag),
        logical(1L))
}
ok <- sim$status == "ok"

## ---- 1. timing ---------------------------------------------------------
timing <- summarize_timing(sim)
timing <- merge(timing, unique(sim[c("scenario_id", "bin_size", "purposes")]),
    by = c("scenario_id", "bin_size"))
timing <- timing[order(timing$n_samples * timing$n_bins), ]
utils::write.csv(timing, file.path(out_dir, "timing.csv"), row.names = FALSE)

## M6: share of the call spent outside lmer + the KR anova (emmeans steps).
time_fit_only <- function(dat) {
    dat <- dat[is.finite(dat$log_norm_cov), ]
    dat$distance_scaled <- dat$distance_from_anchor /
        (max(dat$distance_from_anchor) / 10)
    t0 <- proc.time()[["elapsed"]]
    fit <- suppressMessages(suppressWarnings(lmerTest::lmer(
        log_norm_cov ~ distance_scaled * group + (1 + distance_scaled | sample),
        data = dat, REML = TRUE,
        control = lme4::lmerControl(check.conv.singular = "ignore",
            check.conv.grad = "ignore", check.conv.hess = "ignore"))))
    suppressMessages(suppressWarnings(stats::anova(fit, ddf = "Kenward-Roger")))
    proc.time()[["elapsed"]] - t0
}
m6_cells <- grid[has_tag(grid, "timing_small") | has_tag(grid, "timing_large") |
    has_tag(grid, "timing_large_1kb") |
    (grid$G == 2 & grid$n == 5 & grid$anchor_frac == 0.5 & grid$Delta == 0 &
        grid$tau == 2e-6 & grid$s_ref == 6e-6 & grid$L == 238080), ]
m6 <- do.call(rbind, lapply(seq_len(nrow(m6_cells)), function(i) {
    sc <- scenario_from_row(m6_cells[i, ])
    reps <- seq_len(min(10L, m6_cells$R[i]))
    full <- sim$wall_time[sim$scenario_id == sc$scenario_id &
        sim$bin_size == sc$bin_size & sim$replicate_id %in% reps]
    fit_only <- vapply(reps, function(r) time_fit_only(generate_replicate(sc, r)),
        numeric(1L))
    tibble::tibble(scenario_id = sc$scenario_id, bin_size = sc$bin_size,
        G = sc$G, n = sc$n[1L], L = sc$L, n_fits = length(reps),
        full_median = stats::median(full),
        fit_only_median = stats::median(fit_only),
        emmeans_share = 1 - stats::median(fit_only) / stats::median(full))
}))
utils::write.csv(m6, file.path(out_dir, "emmeans-share.csv"), row.names = FALSE)

## ---- 2. singular rate (M1) and 3. failure rate -------------------------
null_rows <- sim[has_tag(sim, "null") & sim$Delta == 0, ]
sing_by_cell <- do.call(rbind, lapply(split(null_rows,
    paste(null_rows$G, null_rows$n)), function(d) {
    dk <- d[d$status == "ok", ]
    tibble::tibble(G = d$G[1L], n = d$n[1L], n_fits = nrow(d),
        n_ok = nrow(dk), singular_rate = mean(dk$singular),
        corr_at_minus1 = mean(dk$re_corr < -0.999, na.rm = TRUE),
        median_re_corr = stats::median(dk$re_corr, na.rm = TRUE),
        median_df2 = stats::median(dk$df2), n_samples = d$n_samples[1L])
}))
sing_by_cell <- sing_by_cell[order(sing_by_cell$G, sing_by_cell$n), ]
sing_all <- null_rows[null_rows$status == "ok", ]
singular_overall <- mean(sing_all$singular)
singular_ci <- cp_interval(sum(sing_all$singular), nrow(sing_all))
utils::write.csv(sing_by_cell, file.path(out_dir, "singular-by-cell.csv"),
    row.names = FALSE)

fail_by_cell <- do.call(rbind, lapply(split(sim,
    paste(sim$scenario_id, sim$bin_size)), function(d) {
    tibble::tibble(scenario_id = d$scenario_id[1L], bin_size = d$bin_size[1L],
        purposes = d$purposes[1L], G = d$G[1L], n = d$n[1L], L = d$L[1L],
        anchor_frac = NA_real_, n_fits = nrow(d),
        n_failed = sum(d$status != "ok"),
        failure_rate = mean(d$status != "ok"),
        n_conv_nonzero = sum(d$conv_code != 0, na.rm = TRUE),
        n_with_warnings = sum(d$n_warnings > 0, na.rm = TRUE))
}))
fail_by_cell$anchor_frac <- grid$anchor_frac[match(fail_by_cell$scenario_id,
    grid$scenario_id)]
utils::write.csv(fail_by_cell, file.path(out_dir, "failures-by-cell.csv"),
    row.names = FALSE)
fail_msgs <- table(sim$status[sim$status != "ok"])

## ---- 4. M3: Type I at three anchors ------------------------------------
m3_rows <- sim[has_tag(sim, "m3"), ]
m3 <- do.call(rbind, lapply(split(m3_rows, m3_rows$scenario_id), function(d) {
    rej <- sum(d$p_value <= 0.05, na.rm = TRUE)  # failure = non-rejection
    ci <- cp_interval(rej, nrow(d))
    tibble::tibble(scenario_id = d$scenario_id[1L],
        anchor_frac = grid$anchor_frac[match(d$scenario_id[1L],
            grid$scenario_id)],
        n_fits = nrow(d), n_failed = sum(d$status != "ok"),
        n_reject = rej, typeI = rej / nrow(d),
        ci_lower = ci[["lower"]], ci_upper = ci[["upper"]],
        singular_rate = mean(d$singular, na.rm = TRUE))
}))
m3 <- m3[order(m3$anchor_frac), ]
m3_tab <- as.matrix(cbind(m3$n_reject, m3$n_fits - m3$n_reject))
m3_fisher <- stats::fisher.test(m3_tab)$p.value
utils::write.csv(m3, file.path(out_dir, "m3-typeI-by-anchor.csv"),
    row.names = FALSE)

## ---- 5. collapse hypothesis (2.5) ---------------------------------------
col_rows <- sim[grepl("collapse", sim$purposes), ]
col_rows$variant <- ifelse(has_tag(col_rows, "collapse_ref"), "ref",
    ifelse(has_tag(col_rows, "collapse_sref"), "s_ref x2",
        "tau x0.5"))
col_cells <- do.call(rbind, lapply(split(col_rows, col_rows$scenario_id),
    function(d) {
    rej <- sum(d$p_value <= 0.05, na.rm = TRUE)
    ci <- cp_interval(rej, nrow(d))
    tibble::tibble(scenario_id = d$scenario_id[1L], variant = d$variant[1L],
        s_ref = d$s_ref[1L], tau = d$tau[1L], Delta = d$Delta[1L],
        tau_eff = d$tau_eff[1L], delta_std = d$delta_std[1L],
        n_fits = nrow(d), power = rej / nrow(d), ci_lower = ci[["lower"]],
        ci_upper = ci[["upper"]], n_failed = sum(d$status != "ok"),
        n_reject = rej)
}))
col_cells <- col_cells[order(col_cells$delta_std, match(col_cells$variant,
    c("ref", "s_ref x2", "tau x0.5"))), ]
utils::write.csv(col_cells, file.path(out_dir, "collapse-power.csv"),
    row.names = FALSE)

compare_to_ref <- function(variant) {
    ref <- col_cells[col_cells$variant == "ref", ]
    oth <- col_cells[col_cells$variant == variant, ]
    ref <- ref[order(ref$delta_std), ]
    oth <- oth[order(oth$delta_std), ]
    per <- do.call(rbind, lapply(seq_len(nrow(ref)), function(k) {
        pt <- stats::prop.test(c(ref$n_reject[k], oth$n_reject[k]),
            c(ref$n_fits[k], oth$n_fits[k]), correct = FALSE)
        tibble::tibble(variant = variant, delta_std_ref = ref$delta_std[k],
            delta_std_variant = oth$delta_std[k],
            power_ref = ref$power[k], power_variant = oth$power[k],
            diff = ref$power[k] - oth$power[k],
            ci_lower = pt$conf.int[1L], ci_upper = pt$conf.int[2L],
            ci_contains_0 = pt$conf.int[1L] <= 0 && pt$conf.int[2L] >= 0)
    }))
    pt <- stats::prop.test(c(sum(ref$n_reject), sum(oth$n_reject)),
        c(sum(ref$n_fits), sum(oth$n_fits)), correct = FALSE)
    pooled <- tibble::tibble(variant = variant, delta_std_ref = NA_real_,
        delta_std_variant = NA_real_,
        power_ref = sum(ref$n_reject) / sum(ref$n_fits),
        power_variant = sum(oth$n_reject) / sum(oth$n_fits),
        diff = sum(ref$n_reject) / sum(ref$n_fits) -
            sum(oth$n_reject) / sum(oth$n_fits),
        ci_lower = pt$conf.int[1L], ci_upper = pt$conf.int[2L],
        ci_contains_0 = pt$conf.int[1L] <= 0 && pt$conf.int[2L] >= 0)
    rbind(per, pooled)
}
collapse_cmp <- rbind(compare_to_ref("s_ref x2"), compare_to_ref("tau x0.5"))
utils::write.csv(collapse_cmp, file.path(out_dir, "collapse-comparison.csv"),
    row.names = FALSE)

## ---- M4: arithmetic behind the Kruskal-Wallis claim ---------------------
## Perfect separation H and chi-square p, and the exact null probability of
## perfect separation, by enumerating all splits of the ranks.
kw_perfect <- function(n, G) {
    x <- seq_len(n * G)
    g <- rep(seq_len(G), each = n)
    h <- unname(stats::kruskal.test(x, g)$statistic)
    p <- stats::kruskal.test(x, g)$p.value
    if (G == 2L) {
        prob <- 2 / choose(n * G, n)
    } else {
        prob <- NA_real_
    }
    tibble::tibble(G = G, n = n, H_max = h, chisq_p_at_max = p,
        reject_at_0.05 = p < 0.05, null_prob_perfect_separation = prob)
}
m4 <- rbind(kw_perfect(2L, 2L), kw_perfect(3L, 2L), kw_perfect(2L, 3L),
    kw_perfect(3L, 3L))
utils::write.csv(m4, file.path(out_dir, "m4-kw-arithmetic.csv"),
    row.names = FALSE)

## ---- provenance ---------------------------------------------------------
sha <- tryCatch(system2("git", c("-C", "..", "rev-parse", "--short", "HEAD"),
    stdout = TRUE, stderr = FALSE), error = function(e) "unknown")
writeLines(c(paste("package git SHA:", sha), paste("date:", Sys.time()),
    "", utils::capture.output(utils::sessionInfo())),
    file.path(out_dir, "pilot-session-info.txt"))

## ---- report -------------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
md_table <- function(df) {
    cols <- names(df)
    body <- apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "),
        " |"))
    c(paste0("| ", paste(cols, collapse = " | "), " |"),
        paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|"),
        body)
}

tm <- timing
tm_tab <- data.frame(
    design = sprintf("G=%d, n=%d, L=%s kb, %s kb bins", tm$G, tm$n,
        format(tm$L / 1000), format(tm$bin_size / 1000)),
    samples = tm$n_samples, bins_per_sample = tm$n_bins / tm$n_samples,
    fits = tm$n_fits, median_s = fmt(tm$time_median, 2),
    p95_s = fmt(tm$time_p95, 2), max_s = fmt(tm$time_max, 2))
extreme <- tm[has_tag(tm, "timing_small") | has_tag(tm, "timing_large") |
    has_tag(tm, "timing_large_1kb"), ]
extreme_tab <- data.frame(
    design = sprintf("%s: G=%d, n=%d, L=%s kb, %s kb bins",
        sub(";.*", "", extreme$purposes), extreme$G, extreme$n,
        format(extreme$L / 1000), format(extreme$bin_size / 1000)),
    obs = extreme$n_bins, median_s = fmt(extreme$time_median, 2),
    p95_s = fmt(extreme$time_p95, 2))
m6_tab <- data.frame(
    design = sprintf("G=%d, n=%d, L=%s kb, %s kb bins", m6$G, m6$n,
        format(m6$L / 1000), format(m6$bin_size / 1000)),
    full_call_median_s = fmt(m6$full_median, 2),
    fit_and_KR_anova_median_s = fmt(m6$fit_only_median, 2),
    emmeans_share = sprintf("%.0f%%", 100 * m6$emmeans_share))

sing_tab <- data.frame(G = sing_by_cell$G, n = sing_by_cell$n,
    samples = sing_by_cell$n_samples, fits = sing_by_cell$n_fits,
    singular = sprintf("%.0f%%", 100 * sing_by_cell$singular_rate),
    `corr_at_-1` = sprintf("%.0f%%", 100 * sing_by_cell$corr_at_minus1),
    median_df2 = fmt(sing_by_cell$median_df2, 1), check.names = FALSE)
fail_tab <- data.frame(
    cell = sprintf("G=%d n=%d L=%s kb anchor=%s bin=%s kb", fail_by_cell$G,
        fail_by_cell$n, format(fail_by_cell$L / 1000),
        format(fail_by_cell$anchor_frac), format(fail_by_cell$bin_size / 1000)),
    fits = fail_by_cell$n_fits, failed = fail_by_cell$n_failed,
    conv_nonzero = fail_by_cell$n_conv_nonzero)
m3_tab_md <- data.frame(anchor = sprintf("%.1f L", m3$anchor_frac),
    fits = m3$n_fits, rejected = m3$n_reject,
    typeI = fmt(m3$typeI, 2),
    CP95 = sprintf("%.2f to %.2f", m3$ci_lower, m3$ci_upper),
    singular = sprintf("%.0f%%", 100 * m3$singular_rate))
col_md <- data.frame(variant = col_cells$variant,
    s_ref = format(col_cells$s_ref), tau = format(col_cells$tau),
    Delta = formatC(col_cells$Delta, format = "e", digits = 2),
    tau_eff = formatC(col_cells$tau_eff, format = "e", digits = 3),
    `Delta/tau_eff` = fmt(col_cells$delta_std, 3), fits = col_cells$n_fits,
    power = fmt(col_cells$power, 2),
    CP95 = sprintf("%.2f to %.2f", col_cells$ci_lower, col_cells$ci_upper),
    check.names = FALSE)
cmp_md <- data.frame(variant = collapse_cmp$variant,
    `Delta/tau_eff` = ifelse(is.na(collapse_cmp$delta_std_ref), "pooled",
        fmt(collapse_cmp$delta_std_ref, 2)),
    power_ref = fmt(collapse_cmp$power_ref, 2),
    power_variant = fmt(collapse_cmp$power_variant, 2),
    diff = sprintf("%+.2f", collapse_cmp$diff),
    CI95 = sprintf("%+.2f to %+.2f", collapse_cmp$ci_lower,
        collapse_cmp$ci_upper),
    `CI contains 0` = collapse_cmp$ci_contains_0, check.names = FALSE)
m4_md <- data.frame(G = m4$G, n = m4$n, H_max = fmt(m4$H_max, 3),
    chisq_p = fmt(m4$chisq_p_at_max, 4), rejects = m4$reject_at_0.05,
    null_prob_separation = ifelse(is.na(m4$null_prob_perfect_separation),
        "-", fmt(m4$null_prob_perfect_separation, 3)))

report <- c(
    "# Tier 0 pilot report",
    "",
    sprintf("Generated by `scripts/pilot_timing.R` on %s; package git SHA %s.",
        format(Sys.time(), "%Y-%m-%d %H:%M"), sha),
    sprintf(paste("R = %d per cell (M3 cells R = %d, collapse cells R = %d).",
        "%d fits in %.0f s of serial compute. Diagnostic only: at R = 20",
        "every proportion below has a Monte Carlo SE of up to 0.11."),
        R_default, R_m3, R_collapse, nrow(sim),
        sum(sim$wall_time, na.rm = TRUE)),
    "",
    "## 1. Wall time per `test_decay_slopes()` call",
    "",
    "Smallest and largest designs in the grid (the pilot also runs the",
    "largest design at 1 kb bins, the package default):",
    "",
    md_table(extreme_tab),
    "",
    "All timed cells:",
    "",
    md_table(tm_tab),
    "",
    "Share of the call spent outside `lmer` plus the Kenward-Roger anova",
    "(the `emmeans` slope and Tukey steps, M6), first 10 replicates per design:",
    "",
    md_table(m6_tab),
    "",
    "## 2. Singular rate under the baseline null (M1)",
    "",
    sprintf(paste("Overall: %.1f%% of %d successful fits have",
        "`isSingular(fit)` TRUE (Clopper-Pearson 95%% CI %.1f%% to %.1f%%).",
        "Reference geometry, `Delta = 0`, 5 kb bins."),
        100 * singular_overall, nrow(sing_all), 100 * singular_ci[["lower"]],
        100 * singular_ci[["upper"]]),
    sprintf("Median intercept-slope random-effect correlation: %s; fits with correlation at -1 (< -0.999): %.0f%%.",
        fmt(stats::median(sing_all$re_corr, na.rm = TRUE), 3),
        100 * mean(sing_all$re_corr < -0.999, na.rm = TRUE)),
    "",
    md_table(sing_tab),
    "",
    "`median_df2` is the Kenward-Roger denominator df (M2).",
    "",
    "## 3. Failure rate by scenario",
    "",
    sprintf("Overall: %d of %d fits failed (%.1f%%).", sum(!ok), nrow(sim),
        100 * mean(!ok)),
    if (length(fail_msgs)) c("", "Failure messages:", "",
        paste0("- ", names(fail_msgs), " (", as.integer(fail_msgs), ")")),
    "",
    md_table(fail_tab),
    "",
    "## 4. M3 preliminary: Type I at three anchor positions",
    "",
    "Baseline null, `G = 2`, `n = 5`, reference geometry; failure counted as",
    "non-rejection (Q6).",
    "",
    md_table(m3_tab_md),
    "",
    sprintf(paste("Fisher exact test of equal rejection rates across the",
        "three anchors: p = %.3f. At this R the test cannot distinguish",
        "rates that differ by less than about 0.25; a non-significant",
        "result is not evidence of equality."), m3_fisher),
    "",
    "## 5. Collapse hypothesis sanity (PLAN.md 2.5)",
    "",
    "`G = 2`, `n = 5`, two `Delta` levels. In this generator `tau_eff =",
    "sqrt(tau^2 + sigma_bin^2 / Sxx)` does not depend on `s_ref` (a shared",
    "`s_ref` shifts every sample's slope equally, and the median",
    "normalization removes the level), so the `s_ref` pair has a `tau_eff`",
    "ratio of 1 and tests that power does not depend on `s_ref`. The `tau`",
    "pair changes `tau_eff` by a known factor and rescales `Delta` to hold",
    "`Delta / tau_eff` matched; it is the informative test of the collapse.",
    "",
    md_table(col_md),
    "",
    "Difference in power (reference minus variant) at matched `Delta / tau_eff`:",
    "",
    md_table(cmp_md),
    "",
    "## Arithmetic behind M4 (Kruskal-Wallis at small n)",
    "",
    "Perfect separation of the ranks, which gives the largest possible H:",
    "",
    md_table(m4_md),
    "",
    "## Files",
    "",
    "`pilot-replicates.csv` (every fit), `timing.csv`, `emmeans-share.csv`,",
    "`singular-by-cell.csv`, `failures-by-cell.csv`,",
    "`m3-typeI-by-anchor.csv`, `collapse-power.csv`,",
    "`collapse-comparison.csv`, `m4-kw-arithmetic.csv`,",
    "`pilot-session-info.txt`."
)
writeLines(report, file.path(out_dir, "pilot-report.md"))
cat("Wrote", file.path(out_dir, "pilot-report.md"), "\n")
