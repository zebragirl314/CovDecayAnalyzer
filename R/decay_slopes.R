## =========================================================================
## Distance from an anchor position and per-sample decay slope diagnostics.
##
## The formal test of group differences in coverage decay is the linear
## mixed-effects model in decay_tests.R. The per-sample slopes computed
## here are a diagnostic view and the input to plot_decay_slopes().
## =========================================================================

#' Add distance-from-anchor to combined coverage data
#'
#' Computes the absolute distance in base pairs between each bin midpoint
#' and a user-supplied anchor position, and returns the input with a new
#' `distance_from_anchor` column.
#'
#' @details
#'
#' The anchor position is the reference coordinate against which decay
#' will be measured. In HCMV work, this is typically the OriLyt
#' coordinate. To analyze decay from the genome start or end, use
#' `anchor_position = 0` or `anchor_position = genome_length`
#' respectively.
#'
#' Distance is computed as `abs(bin_mid - anchor_position)`. This
#' treats the genome as linear with respect to the anchor. Herpesvirus
#' genomes are linear during the packaged state, so this convention is
#' appropriate. For circular genomes, users would need to compute
#' distance as the shorter of the two arc lengths; that variant is not
#' currently implemented in the package.
#'
#' The returned distances are always non-negative. Sign information
#' (whether a bin is upstream or downstream of the anchor) is discarded,
#' since decay is generally assumed to be symmetric about the anchor.
#' If asymmetric decay is of interest, add a signed distance column
#' manually before fitting.
#'
#' @param combined A data.frame or tibble, typically the output of
#'   [load_all_samples()]. Must contain a numeric `bin_mid` column
#'   giving the midpoint coordinate of each bin.
#' @param anchor_position Non-negative numeric scalar. Reference
#'   coordinate in base pairs from which distances are measured.
#'   Common choices include the OriLyt coordinate (approximately
#'   `93000-94000` for HCMV TB40/E; use the value that matches your
#'   reference annotation), genome start (`0`), or genome end
#'   (`genome_length`).
#' @param genome_length Optional positive integer. If supplied, used to
#'   sanity-check that `anchor_position` does not exceed the reference
#'   length. Recommended in scripts to catch coordinate errors early.
#'   Ignored if `NULL`.
#'
#' @return A tibble with the same columns and rows as `combined`, plus:
#'   \describe{
#'     \item{`distance_from_anchor`}{Numeric. `abs(bin_mid -
#'       anchor_position)`, in base pairs.}
#'   }
#'
#' @seealso [fit_decay_slopes()], [test_decay_slopes()].
#'
#' @examples
#' combined <- tibble::tibble(
#'     sample  = "s1", group = "g1",
#'     bin_mid = seq(500, 4500, by = 1000),
#'     log_norm_cov = rnorm(5)
#' )
#' add_distance_from_anchor(combined,
#'     anchor_position = 2500, genome_length = 5000)
#'
#' @importFrom tibble as_tibble
#' @export
add_distance_from_anchor <- function(combined,
                                    anchor_position,
                                    genome_length = NULL) {
    if (!is.data.frame(combined)) {
        stop("`combined` must be a data.frame or tibble.")
    }
    if (!"bin_mid" %in% names(combined)) {
        stop("`combined` must contain a `bin_mid` column.")
    }
    .validate_position(anchor_position, "anchor_position",
        genome_length = genome_length)

    out <- combined
    out$distance_from_anchor <- abs(out$bin_mid - anchor_position)
    tibble::as_tibble(out)
}


#' Fit per-sample log-linear decay slopes (diagnostic)
#'
#' Fits a simple linear regression per sample of the form
#' `log_norm_cov ~ distance_from_anchor` and returns each sample's slope
#' estimate, standard error, and residual degrees of freedom. This is a
#' diagnostic and visualization aid, not a formal test.
#'
#' @details
#'
#' Per-sample slope estimates are useful for:
#'
#' - Inspecting the spread of decay rates within each treatment group,
#'   e.g. as a boxplot ([plot_decay_slopes()]).
#' - Detecting outlier samples whose decay is qualitatively different
#'   from others in the same group.
#' - Sanity-checking that the decay direction (typically negative slope
#'   in log space) matches expectations.
#'
#' Per-sample slopes are *not* recommended as the basis for a formal
#' group comparison. The two-stage procedure (per-sample slope, then
#' compare across groups) throws away the uncertainty in each per-sample
#' fit, weighting well-fit and poorly-fit samples equally. Use
#' [test_decay_slopes()] instead, which fits a single mixed model to
#' all data at once and preserves within-sample precision information.
#'
#' Samples with fewer than `min_bins` bins after distance filtering are
#' skipped with a warning. This can happen when `window` is small
#' relative to `bin_size` (few bins remain within the window) or when a
#' sample's coverage file yielded few usable bins after normalization.
#'
#' @param combined A data.frame or tibble from [load_all_samples()].
#'   Must contain columns `sample`, `group`, `bin_mid`, and
#'   `log_norm_cov`. If a `distance_from_anchor` column is already
#'   present, `anchor_position` is ignored.
#' @param anchor_position Non-negative numeric scalar. Reference
#'   coordinate for distance calculation. Ignored if `combined`
#'   already has a `distance_from_anchor` column.
#' @param window Optional positive numeric scalar. If supplied, only
#'   bins within `window` base pairs of the anchor (inclusive) are
#'   used in the per-sample fit. Useful for restricting analysis to a
#'   local region around the anchor (e.g. `window = 50000` for a
#'   50 kb window around OriLyt). If `NULL`, all bins in `combined`
#'   are used.
#' @param min_bins Positive integer. Minimum number of usable bins
#'   required to attempt a fit for a given sample. Samples with fewer
#'   than this are skipped with a warning. Default 5.
#'
#' @return A tibble with one row per successfully-fit sample:
#'   \describe{
#'     \item{`sample`}{Character. Sample name.}
#'     \item{`group`}{Character. Group label.}
#'     \item{`slope`}{Numeric. Slope estimate for
#'       `distance_from_anchor` (units: change in `log_norm_cov` per
#'       base pair). Typically negative when coverage decays away from
#'       the anchor.}
#'     \item{`se_slope`}{Numeric. Standard error of the slope estimate.}
#'     \item{`df`}{Integer. Residual degrees of freedom of the sample's
#'       linear fit (approximately `n_bins - 2`).}
#'     \item{`n_bins`}{Integer. Number of bins used in the fit after
#'       applying `window` and dropping non-finite rows.}
#'   }
#'
#' @seealso [test_decay_slopes()] for the recommended formal group
#'   comparison; [plot_decay_slopes()] for a per-sample slope
#'   visualization.
#'
#' @examples
#' set.seed(1)
#' fake <- do.call(rbind, lapply(paste0("s", 1:4), function(sname) {
#'     tibble::tibble(
#'         sample = sname,
#'         group  = if (sname %in% c("s1", "s2")) "g1" else "g2",
#'         bin_mid = seq(500, 100000, by = 1000),
#'         log_norm_cov = -1e-5 * seq(500, 100000, by = 1000) +
#'             rnorm(100, sd = 0.1)
#'     )
#' }))
#' fit_decay_slopes(fake, anchor_position = 0, window = 50000)
#'
#' @importFrom stats lm coef vcov
#' @importFrom tibble tibble
#' @export
fit_decay_slopes <- function(combined,
                            anchor_position = NULL,
                            window = NULL,
                            min_bins = 5L) {
    if (!is.data.frame(combined)) {
        stop("`combined` must be a data.frame or tibble.")
    }
    required <- c("sample", "group", "bin_mid", "log_norm_cov")
    missing_cols <- setdiff(required, names(combined))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "`combined` is missing required column(s): %s",
            paste(missing_cols, collapse = ", ")
        ))
    }
    .validate_positive_integer(min_bins, "min_bins")

    if (!"distance_from_anchor" %in% names(combined)) {
        if (is.null(anchor_position)) {
            stop("Provide `anchor_position` or a `distance_from_anchor` ",
                "column in `combined`.")
        }
        combined <- add_distance_from_anchor(combined, anchor_position)
    }

    if (!is.null(window)) {
        if (!is.numeric(window) || length(window) != 1L ||
            !is.finite(window) || window <= 0) {
            stop("`window` must be a single positive finite number.")
        }
        combined <- combined[combined$distance_from_anchor <= window, ,
            drop = FALSE]
    }

    ## Drop bins where the response is not finite (e.g. samples that
    ## failed normalization due to zero median coverage).
    combined <- combined[is.finite(combined$log_norm_cov) &
        is.finite(combined$distance_from_anchor), , drop = FALSE]

    samples <- unique(combined$sample)
    if (length(samples) == 0L) {
        stop("No samples remain after filtering. ",
            "Check `window` and input data.")
    }

    results <- lapply(samples, function(sname) {
        sdat <- combined[combined$sample == sname, , drop = FALSE]
        if (nrow(sdat) < min_bins) {
            warning(sprintf(
                "Sample '%s' has only %d bins after filtering; skipped.",
                sname, nrow(sdat)
            ))
            return(NULL)
        }
        fit <- tryCatch(
            stats::lm(log_norm_cov ~ distance_from_anchor, data = sdat),
            error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)

        cf <- stats::coef(fit)
        se <- sqrt(diag(stats::vcov(fit)))
        tibble::tibble(
            sample   = sname,
            group    = sdat$group[1],
            slope    = unname(cf["distance_from_anchor"]),
            se_slope = unname(se["distance_from_anchor"]),
            df       = fit$df.residual,
            n_bins   = nrow(sdat)
        )
    })

    results <- results[!vapply(results, is.null, logical(1L))]
    if (length(results) == 0L) {
        stop("No sample yielded a valid slope fit.")
    }
    do.call(rbind, results)
}
