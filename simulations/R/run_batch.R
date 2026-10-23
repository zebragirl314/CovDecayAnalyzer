## =========================================================================
## Batch runner: one (scenario, bin size) x a block of replicates.
##
## Minimal in A2-1 (primary LMM arm only). A2-2 adds the benchmark arms to
## this loop so every arm is paired on the same dataset.
## =========================================================================

#' Run a batch of replicates.
#'
#' @param batch One row of [expand_batches()]: the scenario columns plus
#'   `rep_start` and `rep_end`.
#' @param methods_files Character vector of the package files that define
#'   the statistics. Unused in the body; it exists so that a change to those
#'   files invalidates the target (PLAN.md 5.1).
#' @return A tibble, one row per replicate, with scenario parameters
#'   (`delta_std`, `tau_eff`, ...) and the [analyze_lmm()] columns.
run_batch <- function(batch, methods_files = NULL) {
    ensure_package()
    sc <- scenario_from_row(batch)
    reps <- seq.int(batch$rep_start, batch$rep_end)
    rows <- lapply(reps, function(r) {
        dat <- generate_replicate(sc, r)
        analyze_lmm(dat)
    })
    res <- do.call(rbind, rows)
    meta <- scenario_to_row(sc)
    meta$scenario_id <- NULL
    meta$bin_size <- NULL
    cbind(res[c("scenario_id", "replicate_id", "bin_size")], meta,
        res[setdiff(names(res), c("scenario_id", "replicate_id", "bin_size"))])
}
