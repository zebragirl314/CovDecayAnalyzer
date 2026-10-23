## =========================================================================
## Primary analysis arm: wrapper around `test_decay_slopes()` (PLAN.md 2.4).
##
## Returns one row per call, never an error: a failed fit comes back with
## `status` set to the error message and NA results, so summaries can count
## failures (unconditional rate, Q6) instead of silently dropping rows.
## =========================================================================

MAX_GROUPS <- 4L

.empty_lmm_row <- function(scenario_id, replicate_id, bin_size) {
    slopes <- stats::setNames(as.list(rep(NA_real_, MAX_GROUPS)),
        paste0("slope_g", seq_len(MAX_GROUPS)))
    tibble::tibble(
        scenario_id = scenario_id, replicate_id = replicate_id,
        bin_size = bin_size, arm = "lmm",
        status = NA_character_,
        p_value = NA_real_, F_stat = NA_real_, df1 = NA_real_,
        df2 = NA_real_, singular = NA, conv_code = NA_integer_,
        sigma = NA_real_, sd_slope = NA_real_, sd_intercept = NA_real_,
        re_corr = NA_real_,
        tibble::as_tibble(slopes),
        n_samples = NA_integer_, n_bins = NA_integer_,
        n_warnings = NA_integer_, wall_time = NA_real_
    )
}

#' Run the primary LMM test on one simulated dataset.
#'
#' @param data Tibble from [generate_replicate()].
#' @param scenario_id,replicate_id,bin_size Identifiers copied into the row;
#'   default to the attributes of `data`.
#' @return One-row tibble: `scenario_id`, `replicate_id`, `bin_size`, `arm`,
#'   `status` (`"ok"` or the error message), `p_value`, `F_stat`, `df1`,
#'   `df2`, `singular`, `conv_code`, `sigma` (residual SD), `sd_slope` and
#'   `sd_intercept` (random-effect SDs; `sd_slope` per bp), `re_corr`
#'   (intercept-slope correlation), `slope_g1` to `slope_g4` (estimated
#'   group slopes per bp on the LMM scale, i.e. `-s`; NA beyond `G`),
#'   `n_samples`, `n_bins`, `n_warnings`, and `wall_time` (seconds for the
#'   `test_decay_slopes()` call, including its `emmeans` steps).
analyze_lmm <- function(data,
                        scenario_id = attr(data, "scenario_id"),
                        replicate_id = attr(data, "replicate_id"),
                        bin_size = attr(data, "bin_size")) {
    ensure_package()
    row <- .empty_lmm_row(scenario_id %||% NA_character_,
        replicate_id %||% NA_integer_, bin_size %||% NA_integer_)
    n_warn <- 0L
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(
        withCallingHandlers(
            CovDecayAnalyzer::test_decay_slopes(data),
            warning = function(w) {
                n_warn <<- n_warn + 1L
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) e
    )
    row$wall_time <- proc.time()[["elapsed"]] - t0
    row$n_warnings <- n_warn

    if (inherits(res, "error")) {
        row$status <- conditionMessage(res)
        return(row)
    }

    ext <- tryCatch({
        fit <- res$model
        vc <- lme4::VarCorr(fit)$sample
        sds <- attr(vc, "stddev")
        cr <- attr(vc, "correlation")
        slope_col <- grep("\\.trend$", names(res$group_slopes), value = TRUE)
        sl <- res$group_slopes[[slope_col[1L]]]
        list(
            sigma = stats::sigma(fit),
            sd_slope = unname(sds[length(sds)]) / res$distance_scale,
            sd_intercept = unname(sds[1L]),
            re_corr = if (length(sds) == 2L) unname(cr[1L, 2L]) else NA_real_,
            slopes = sl
        )
    }, error = function(e) e)

    row$p_value <- res$p_value
    row$F_stat <- res$F_stat
    row$df1 <- res$df1
    row$df2 <- res$df2
    row$singular <- res$singular
    row$conv_code <- if (is.null(res$conv_code)) NA_integer_ else res$conv_code
    row$n_samples <- as.integer(res$n_samples)
    row$n_bins <- as.integer(res$n_bins)

    if (inherits(ext, "error")) {
        row$status <- paste("extract:", conditionMessage(ext))
        return(row)
    }
    row$sigma <- ext$sigma
    row$sd_slope <- ext$sd_slope
    row$sd_intercept <- ext$sd_intercept
    row$re_corr <- ext$re_corr
    k <- min(length(ext$slopes), MAX_GROUPS)
    for (g in seq_len(k)) {
        row[[paste0("slope_g", g)]] <- ext$slopes[g]
    }
    row$status <- "ok"
    row
}

`%||%` <- function(x, y) if (is.null(x)) y else x
