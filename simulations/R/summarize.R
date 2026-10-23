## =========================================================================
## Summaries. A2-1 defines only what Tier 0 needs; the Type I and power
## summaries (Clopper-Pearson, Q7) are completed in A2-2.
## =========================================================================

## Exact (Clopper-Pearson) interval for a proportion.
cp_interval <- function(x, n, level = 0.95) {
    if (n == 0L) return(c(lower = NA_real_, upper = NA_real_))
    ci <- stats::binom.test(x, n, conf.level = level)$conf.int
    c(lower = ci[1L], upper = ci[2L])
}

#' Wall-time summary by scenario and bin size.
#'
#' @param sim Tibble of replicate rows from [run_batch()].
#' @return Tibble: one row per scenario and bin size with the number of
#'   fits, the number of bins and samples, and the median, 95th percentile,
#'   and maximum wall time in seconds.
summarize_timing <- function(sim) {
    key <- paste(sim$scenario_id, sim$bin_size)
    out <- lapply(split(sim, key), function(d) {
        tibble::tibble(
            scenario_id = d$scenario_id[1L], bin_size = d$bin_size[1L],
            G = d$G[1L], n = d$n[1L], L = d$L[1L],
            n_samples = d$n_samples[1L], n_bins = d$n_bins[1L],
            n_fits = nrow(d),
            time_median = stats::median(d$wall_time, na.rm = TRUE),
            time_p95 = unname(stats::quantile(d$wall_time, 0.95,
                na.rm = TRUE)),
            time_max = max(d$wall_time, na.rm = TRUE)
        )
    })
    do.call(rbind, unname(out))
}
