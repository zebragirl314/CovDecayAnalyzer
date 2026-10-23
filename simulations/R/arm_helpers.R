## =========================================================================
## Shared helpers for the analysis arms added in A2-2 prep. Every arm returns
## rows with the column shape of `analyze_lmm()` (see `.empty_lmm_row()`), so
## results from different arms `rbind` directly. Arms with several outputs
## (KW omnibus, Dunn ungated, Dunn gated) return one row each, distinguished
## by `arm`.
## =========================================================================

.arm_row <- function(data, arm, scenario_id, replicate_id, bin_size) {
    row <- .empty_lmm_row(scenario_id %||% NA_character_,
        replicate_id %||% NA_integer_, bin_size %||% NA_integer_)
    row$arm <- arm
    row
}

## Per-sample OLS slopes via the package's `fit_decay_slopes()`.
.sample_slopes <- function(data) {
    sl <- as.data.frame(CovDecayAnalyzer::fit_decay_slopes(data))
    sl$group <- factor(sl$group)
    sl
}

.group_slope_means <- function(sl) {
    m <- tapply(sl$slope, sl$group, mean)
    out <- rep(NA_real_, MAX_GROUPS)
    out[seq_along(m)] <- as.numeric(m)
    out
}

## Failure rows for every arm name, so a caught error still yields the full
## set of rows that a success would.
.failed_rows <- function(arms, msg, ...) {
    do.call(rbind, lapply(arms, function(a) {
        r <- .arm_row(NULL, a, ...)
        r$status <- msg
        r
    }))
}

## Dunn's test on pooled ranks (normal approximation, no tie correction:
## slopes are continuous), Holm-adjusted. Returns the adjusted p-values.
.dunn_holm <- function(slope, group) {
    N <- length(slope)
    rk <- rank(slope)
    rbar <- tapply(rk, group, mean)
    ng <- tapply(rk, group, length)
    pr <- utils::combn(names(rbar), 2L)
    z <- apply(pr, 2L, function(g) {
        (rbar[[g[1L]]] - rbar[[g[2L]]]) /
            sqrt(N * (N + 1) / 12 * (1 / ng[[g[1L]]] + 1 / ng[[g[2L]]]))
    })
    stats::p.adjust(2 * stats::pnorm(-abs(z)), method = "holm")
}

## Assemble the KW rows: omnibus, Dunn-Holm ungated, Dunn-Holm gated (normal
## approximation), and, when `kw$dunn_perm` is present, the permutation-Dunn
## ungated and gated rows. For the Dunn rows `p_value` is the smallest
## Holm-adjusted pairwise p (the family-wise rejection indicator at any
## alpha); gated versions are 1 unless the KW p-value is below 0.05. Per-pair
## p-values are not stored. `kw_perm_method` and `dunn_perm_method` record
## whether the permutation p-values came from enumeration or random draws.
.kw_rows <- function(prefix, data, kw_fun, scenario_id, replicate_id, bin_size,
                     perm_dunn = FALSE) {
    suffix <- c("", "_dunn_ungated", "_dunn_gated")
    if (perm_dunn) suffix <- c(suffix, "_dunn_perm_ungated", "_dunn_perm_gated")
    arms <- paste0(prefix, suffix)
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch({
        sl <- .sample_slopes(data)
        kw <- kw_fun(sl$slope, sl$group)
        pa <- .dunn_holm(sl$slope, sl$group)
        rows <- lapply(arms, .arm_row, data = data, scenario_id = scenario_id,
            replicate_id = replicate_id, bin_size = bin_size)
        G <- nlevels(sl$group)
        for (i in seq_along(rows)) {
            rows[[i]]$status <- "ok"
            rows[[i]]$n_samples <- nrow(sl)
            rows[[i]]$n_bins <- nrow(data)
            rows[[i]]$n_warnings <- 0L
            rows[[i]]$df1 <- G - 1
            for (g in seq_len(min(G, MAX_GROUPS))) {
                rows[[i]][[paste0("slope_g", g)]] <- .group_slope_means(sl)[g]
            }
        }
        if (!is.null(kw$method)) {
            for (i in seq_along(rows)) {
                rows[[i]]$kw_perm_method <- kw$method
                rows[[i]]$dunn_perm_method <- NA_character_
            }
        }
        rows[[1L]]$p_value <- kw$p
        rows[[1L]]$F_stat <- kw$H
        rows[[2L]]$p_value <- min(pa)
        rows[[3L]]$p_value <- if (kw$p < 0.05) min(pa) else 1
        if (perm_dunn) {
            pp <- min(kw$dunn_perm)
            rows[[4L]]$p_value <- pp
            rows[[5L]]$p_value <- if (kw$p < 0.05) pp else 1
            rows[[4L]]$dunn_perm_method <- kw$method
            rows[[5L]]$dunn_perm_method <- kw$method
        }
        do.call(rbind, rows)
    }, error = function(e) {
        .failed_rows(arms, conditionMessage(e), scenario_id, replicate_id,
            bin_size)
    })
    out$wall_time <- proc.time()[["elapsed"]] - t0
    out
}
