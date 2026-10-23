## =========================================================================
## Batch runner: one (scenario, bin size) x a block of replicates.
##
## Each replicate is generated once and analyzed by all six arms in a single
## pass, so arms are paired by construction (PLAN.md 5.1). The result is a
## long tibble, one row per (replicate, arm, test_variant).
## =========================================================================

## Arm functions in the order they run. The slope-only LMM is a cautionary
## comparator, not a validation target (it ignores real between-sample
## intercept variance and is badly anti-conservative); its output name carries
## the flag so no table can be read as recommending it.
ARM_FUNS <- list(
    lmm = "analyze_lmm",
    lmm_slope_only = "analyze_lmm_slope_only",
    kw_chisq = "analyze_kw_chisq",
    kw_perm = "analyze_kw_perm",
    two_stage = "analyze_two_stage",
    two_stage_pooled = "analyze_two_stage_pooled"
)
ARM_OUTPUT_NAMES <- c(lmm = "lmm", lmm_slope_only = "lmm_slope_only_CAUTIONARY",
    kw_chisq = "kw_chisq", kw_perm = "kw_perm", two_stage = "two_stage",
    two_stage_pooled = "two_stage_pooled")

## Split the internal `arm` string ("kw_chisq_dunn_gated") into the reported
## arm and test_variant ("omnibus", "dunn_ungated", "dunn_gated").
.split_arm <- function(arm) {
    variant <- ifelse(grepl("_dunn_perm_ungated$", arm), "dunn_perm_ungated",
        ifelse(grepl("_dunn_perm_gated$", arm), "dunn_perm_gated",
        ifelse(grepl("_dunn_ungated$", arm), "dunn_ungated",
        ifelse(grepl("_dunn_gated$", arm), "dunn_gated", "omnibus"))))
    base <- sub("_dunn_(perm_)?(un)?gated$", "", arm)
    out <- unname(ARM_OUTPUT_NAMES[base])
    out[is.na(out)] <- base[is.na(out)]
    list(arm = out, test_variant = variant)
}

## Run one arm; never errors. A caught error yields NA result rows with the
## message in `status` (one row per variant the arm would have produced).
.run_arm <- function(name, data, scenario_id, replicate_id, bin_size) {
    fun <- get(ARM_FUN_NAMES[[name]], mode = "function")
    tryCatch(
        fun(data, scenario_id = scenario_id, replicate_id = replicate_id,
            bin_size = bin_size),
        error = function(e) {
            arms <- if (grepl("^kw_", name)) {
                paste0(name, c("", "_dunn_ungated", "_dunn_gated",
                    if (name == "kw_perm") c("_dunn_perm_ungated",
                        "_dunn_perm_gated")))
            } else name
            .failed_rows(arms, conditionMessage(e), scenario_id,
                replicate_id, bin_size)
        }
    )
}
ARM_FUN_NAMES <- unlist(ARM_FUNS)

.compose_long <- function(arm_rows) {
    ## Arms add different optional columns (kw_perm_method); align them.
    cols <- unique(unlist(lapply(arm_rows, names)))
    arm_rows <- lapply(arm_rows, function(r) {
        for (cn in setdiff(cols, names(r))) r[[cn]] <- NA
        r[cols]
    })
    out <- do.call(rbind, arm_rows)
    for (cn in c("kw_perm_method", "dunn_perm_method")) {
        out[[cn]] <- as.character(if (is.null(out[[cn]])) rep(NA, nrow(out))
            else out[[cn]])
    }
    sp <- .split_arm(out$arm)
    out$arm <- sp$arm
    out$test_variant <- sp$test_variant
    out
}

#' Run a batch of replicates through all six analysis arms.
#'
#' @param batch One row of [expand_batches()]: the scenario columns plus
#'   `rep_start` and `rep_end` (and, for Tier 1, `family`, `sweep_factor`,
#'   `sweep_level`).
#' @param methods_files Character vector of the package files that define
#'   the statistics. Unused in the body; it exists so that a change to those
#'   files invalidates the target (PLAN.md 5.1).
#' @return A tibble, one row per (replicate, arm, test_variant), with the
#'   scenario parameters, `arm`, `test_variant`, `kw_perm_method` (`exact_enum`
#'   or `random_perm` on `kw_perm` rows, NA elsewhere), `dunn_perm_method` (same
#'   values, on the `dunn_perm_*` rows only) and the columns of
#'   [analyze_lmm()]. `status` is `"ok"` or the error message of a failed arm,
#'   whose result columns are then NA. `wall_time` on the three rows of one KW
#'   arm is the shared time of that arm.
run_batch <- function(batch, methods_files = NULL) {
    ensure_package()
    sc <- scenario_from_row(batch)
    reps <- seq.int(batch$rep_start, batch$rep_end)
    rows <- lapply(reps, function(r) {
        dat <- tryCatch(generate_replicate(sc, r), error = function(e) e)
        if (inherits(dat, "error")) {
            all_arms <- c("lmm", "lmm_slope_only", "kw_chisq",
                paste0("kw_perm", c("", "_dunn_ungated", "_dunn_gated",
                    "_dunn_perm_ungated", "_dunn_perm_gated")),
                paste0("kw_chisq", c("_dunn_ungated", "_dunn_gated")),
                "two_stage", "two_stage_pooled")
            res <- .failed_rows(all_arms,
                paste("generate:", conditionMessage(dat)), sc$scenario_id, r,
                sc$bin_size)
            return(.compose_long(list(res)))
        }
        .compose_long(lapply(names(ARM_FUNS), .run_arm, data = dat,
            scenario_id = sc$scenario_id, replicate_id = r,
            bin_size = sc$bin_size))
    })
    res <- do.call(rbind, rows)
    meta <- scenario_to_row(sc)
    meta$scenario_id <- NULL
    meta$bin_size <- NULL
    extra <- intersect(c("family", "sweep_factor", "sweep_level"),
        names(batch))
    for (e in extra) meta[[e]] <- as.character(batch[[e]])
    first <- c("scenario_id", "replicate_id", "bin_size", "arm",
        "test_variant")
    cbind(res[first], meta, res[setdiff(names(res), first)])
}
