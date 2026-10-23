## Two-stage parametric arm with a pooled-variance second stage: per-sample
## OLS slope (`fit_decay_slopes`), then a pooled-variance t-test (G = 2) or
## classical one-way ANOVA (G >= 3). Returns one row, arm = "two_stage_pooled".
## For G = 2, `F_stat` is t^2 and `df2` is N - 2. This is the same F test as
## the interaction F of a balanced LMM on per-sample slopes. For G >= 3 it is
## identical to `analyze_two_stage()`, which already uses the classical
## ANOVA; only G = 2 differs (Welch there).
analyze_two_stage_pooled <- function(data,
                                     scenario_id = attr(data, "scenario_id"),
                                     replicate_id = attr(data, "replicate_id"),
                                     bin_size = attr(data, "bin_size")) {
    ensure_package()
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch({
        sl <- .sample_slopes(data)
        G <- nlevels(sl$group)
        row <- .arm_row(data, "two_stage_pooled", scenario_id, replicate_id,
            bin_size)
        if (G == 2L) {
            tt <- stats::t.test(slope ~ group, data = sl, var.equal = TRUE)
            row$F_stat <- unname(tt$statistic)^2
            row$df1 <- 1
            row$df2 <- unname(tt$parameter)
            row$p_value <- tt$p.value
        } else {
            ow <- stats::oneway.test(slope ~ group, data = sl,
                var.equal = TRUE)
            row$F_stat <- unname(ow$statistic)
            row$df1 <- unname(ow$parameter[1L])
            row$df2 <- unname(ow$parameter[2L])
            row$p_value <- ow$p.value
        }
        for (g in seq_len(min(G, MAX_GROUPS))) {
            row[[paste0("slope_g", g)]] <- .group_slope_means(sl)[g]
        }
        row$n_samples <- nrow(sl)
        row$n_bins <- nrow(data)
        row$n_warnings <- 0L
        row$status <- "ok"
        row
    }, error = function(e) {
        .failed_rows("two_stage_pooled", conditionMessage(e), scenario_id,
            replicate_id, bin_size)
    })
    out$wall_time <- proc.time()[["elapsed"]] - t0
    out
}
