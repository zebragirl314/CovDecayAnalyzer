## Random-slope-only LMM comparator (M1): same fixed effects and Kenward-Roger
## F test as the primary model, random effect `(0 + distance_scaled | sample)`
## (no random intercept). Fitted directly with lmerTest rather than through
## `test_decay_slopes()`, so no emmeans step; distance is rescaled the same
## way (max to 10). Returns one row, arm = "lmm_slope_only". `sd_intercept`
## and `re_corr` are NA by construction; slopes are fixed-effect group slopes
## per bp on the LMM scale (-s).
analyze_lmm_slope_only <- function(data,
                                   scenario_id = attr(data, "scenario_id"),
                                   replicate_id = attr(data, "replicate_id"),
                                   bin_size = attr(data, "bin_size")) {
    ensure_package()
    row <- .arm_row(data, "lmm_slope_only", scenario_id, replicate_id,
        bin_size)
    n_warn <- 0L
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(withCallingHandlers({
        d <- data[is.finite(data$log_norm_cov) &
            is.finite(data$distance_from_anchor), ]
        d$group <- factor(d$group)
        scale <- max(d$distance_from_anchor) / 10
        d$distance_scaled <- d$distance_from_anchor / scale
        fit <- suppressMessages(lmerTest::lmer(
            log_norm_cov ~ distance_scaled * group +
                (0 + distance_scaled | sample),
            data = d, REML = TRUE,
            control = lme4::lmerControl(check.conv.singular = "ignore",
                check.conv.grad = "ignore", check.conv.hess = "ignore")))
        an <- stats::anova(fit, ddf = "Kenward-Roger")
        list(fit = fit, an = an, scale = scale, d = d)
    }, warning = function(w) {
        n_warn <<- n_warn + 1L
        invokeRestart("muffleWarning")
    }), error = function(e) e)
    row$wall_time <- proc.time()[["elapsed"]] - t0
    row$n_warnings <- n_warn
    if (inherits(res, "error")) {
        row$status <- conditionMessage(res)
        return(row)
    }
    r <- grep("distance_scaled:group", rownames(res$an))
    row$p_value <- res$an[r, "Pr(>F)"]
    row$F_stat <- res$an[r, "F value"]
    row$df1 <- res$an[r, "NumDF"]
    row$df2 <- res$an[r, "DenDF"]
    row$singular <- lme4::isSingular(res$fit)
    cc <- res$fit@optinfo$conv$opt
    row$conv_code <- if (is.null(cc)) NA_integer_ else as.integer(cc)
    row$sigma <- stats::sigma(res$fit)
    vc <- lme4::VarCorr(res$fit)$sample
    row$sd_slope <- unname(attr(vc, "stddev")[1L]) / res$scale
    b <- lme4::fixef(res$fit)
    G <- nlevels(res$d$group)
    sl <- b[["distance_scaled"]]
    if (G > 1L) {
        sl <- c(sl, sl + b[grep("^distance_scaled:group", names(b))])
    }
    for (g in seq_len(min(G, MAX_GROUPS))) {
        row[[paste0("slope_g", g)]] <- sl[g] / res$scale
    }
    row$n_samples <- length(unique(res$d$sample))
    row$n_bins <- nrow(res$d)
    row$status <- "ok"
    row
}
