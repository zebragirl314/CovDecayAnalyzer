## =========================================================================
## Random-anchor null distribution and empirical p-value.
##
## Refits the linear mixed model at N random genomic anchor positions to
## build a null distribution of Kenward-Roger F-statistics for the
## distance-by-group interaction. The observed F-statistic at the anchor
## of scientific interest (e.g. OriLyt, genome start, genome end) is
## compared to this null to compute an empirical p-value.
##
## This addresses the question "is the group effect at this specific
## anchor stronger than what we would observe at a random genomic
## position?" and controls for the concern that any arbitrary anchor
## choice inflates false positives.
## =========================================================================


#' Build a null distribution of interaction F-statistics from random anchors
#'
#' Repeatedly refits the linear mixed model of [test_decay_slopes()] at
#' random anchor positions drawn uniformly from the genome, and returns
#' the vector of Kenward-Roger F-statistics for the distance-by-group
#' interaction at each random anchor. This vector is the null against
#' which the observed F at a scientifically-motivated anchor (e.g.
#' OriLyt) can be compared.
#'
#' @details
#'
#' The rationale is that a decay-vs-group interaction can arise at any
#' genomic position for many reasons (assembly artifacts, repetitive
#' regions, gene expression differences that happen to correlate with
#' treatment). To attribute a significant interaction specifically to
#' the anchor of interest, we need to show that the effect at that
#' anchor is stronger than at random positions elsewhere in the genome.
#' This is exactly what the null distribution provides.
#'
#' Each iteration:
#'
#' 1. Draws a random position uniformly from `[1, genome_length]`
#'    (optionally excluding a region near the observed anchor).
#' 2. Recomputes `distance_from_anchor` for every bin relative to
#'    that random position.
#' 3. Refits the LMM and extracts the Kenward-Roger F-statistic for the
#'    `distance_from_anchor:group` interaction.
#'
#' Failed fits (singular models, non-converging optimizers) are dropped
#' silently. The returned vector length may therefore be smaller than
#' `n_random`; the actual retention rate is a useful diagnostic and is
#' reported by [test_decay_significance()].
#'
#' Runtime scales linearly in `n_random` times the cost of one LMM fit.
#' On typical HCMV data (~30 samples, ~240 bins each), one fit takes on
#' the order of one to a few seconds, so `n_random = 500` on one core
#' takes roughly 10 to 30 minutes. Parallelize with `BPPARAM =
#' BiocParallel::MulticoreParam(workers = N)`.
#'
#' @param combined A data.frame or tibble from [load_all_samples()].
#'   Must contain columns `sample`, `group`, `bin_mid`, and
#'   `log_norm_cov`.
#' @param genome_length Positive integer. Reference genome length in
#'   base pairs. Random anchor positions are drawn uniformly from
#'   `[1, genome_length]`.
#' @param n_random Positive integer. Number of random anchor positions.
#'   Default 500. The minimum achievable empirical p-value is
#'   approximately `1 / (n_random + 1)`, so higher values give finer
#'   resolution at the tail; 500 is a reasonable default for typical
#'   experimental screens, 1000 or more for publication-grade
#'   significance claims.
#' @param window Optional positive numeric scalar. Passed through
#'   unchanged to each random-anchor call of [test_decay_slopes()].
#'   Use the same value that will be used at the observed anchor of
#'   interest, so the null and observed statistics come from
#'   comparable fits.
#' @param exclude Optional length-2 numeric vector `c(low, high)`. If
#'   supplied, random positions falling in `[low, high]` are rejected.
#'   Useful for excluding a band around a known anchor of interest so
#'   the null does not include overlapping windows. Ignored if `NULL`.
#' @param seed Optional integer. Sets `set.seed()` for reproducibility
#'   of the sampled positions. If `NULL`, no seed is set.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` instance
#'   controlling parallelism. Default [BiocParallel::SerialParam()].
#'   Use [BiocParallel::MulticoreParam] on Unix-like systems or
#'   [BiocParallel::SnowParam] on Windows.
#'
#' @return A numeric vector of length at most `n_random` containing
#'   the Kenward-Roger F-statistics from successful fits at the
#'   sampled random anchor positions. Non-finite values from failed
#'   fits are removed.
#'
#' @references
#' Phipson B, Smyth GK (2010). Permutation p-values should never be
#' zero: calculating exact p-values when permutations are randomly
#' drawn. *Statistical Applications in Genetics and Molecular Biology*
#' 9(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' @seealso [empirical_pvalue()] for the p-value formula that consumes
#'   this null; [test_decay_significance()] for the end-to-end wrapper.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' fake <- do.call(rbind, lapply(paste0("s", 1:8), function(sname) {
#'     grp <- if (as.integer(sub("s", "", sname)) <= 4) "g1" else "g2"
#'     tibble::tibble(
#'         sample = sname, group = grp,
#'         bin_mid = seq(500, 100000, by = 1000),
#'         log_norm_cov = rnorm(100, sd = 0.5)
#'     )
#' }))
#' null_F <- null_distribution_F(
#'     fake, genome_length = 100000, n_random = 10, seed = 1
#' )
#' summary(null_F)
#' }
#'
#' @importFrom BiocParallel bplapply SerialParam
#' @export
null_distribution_F <- function(combined,
                                genome_length,
                                n_random = 500L,
                                window = NULL,
                                exclude = NULL,
                                seed = NULL,
                                BPPARAM = BiocParallel::SerialParam()) {
    if (!is.data.frame(combined)) {
        stop("`combined` must be a data.frame or tibble.")
    }
    .validate_positive_integer(genome_length, "genome_length")
    .validate_positive_integer(n_random, "n_random")

    if (!is.null(exclude)) {
        if (!is.numeric(exclude) || length(exclude) != 2L ||
            any(!is.finite(exclude)) || exclude[1] >= exclude[2]) {
            stop("`exclude` must be a length-2 numeric vector c(low, high).")
        }
    }

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L) {
            stop("`seed` must be a single integer or NULL.")
        }
        set.seed(as.integer(seed))
    }

    ## Sample positions, avoiding the excluded region if requested.
    candidate <- sample.int(genome_length, size = n_random, replace = FALSE)
    if (!is.null(exclude)) {
        candidate <- candidate[
            candidate < exclude[1] | candidate > exclude[2]
        ]
        ## Top up if we removed too many.
        while (length(candidate) < n_random) {
            extra <- sample.int(genome_length,
                size = n_random - length(candidate), replace = FALSE)
            extra <- extra[extra < exclude[1] | extra > exclude[2]]
            candidate <- unique(c(candidate, extra))
            if (length(candidate) >= n_random) break
        }
        candidate <- candidate[seq_len(min(n_random, length(candidate)))]
    }

    F_list <- BiocParallel::bplapply(
        candidate,
        function(pos) {
            tryCatch(
                {
                    res <- test_decay_slopes(
                        combined,
                        anchor_position = pos,
                        window          = window,
                        quiet           = TRUE
                    )
                    res$F_stat
                },
                error = function(e) NA_real_
            )
        },
        BPPARAM = BPPARAM
    )
    F_vec <- unlist(F_list, use.names = FALSE)
    F_vec[is.finite(F_vec)]
}


#' Empirical (Monte Carlo) p-value with continuity correction
#'
#' Computes a right-tailed empirical p-value from an observed test
#' statistic and a null distribution obtained by resampling, using the
#' `(k + 1) / (n + 1)` continuity correction of Phipson & Smyth (2010).
#'
#' @details
#'
#' The formula is
#'
#' \deqn{p = (k + 1) / (n + 1)}{p = (k + 1) / (n + 1)}
#'
#' where `k` is the number of null values greater than or equal to the
#' observed statistic and `n` is the size of the null distribution.
#'
#' The `+1` in both numerator and denominator is not an ad-hoc smoothing
#' term; it is the correct exact p-value expression for Monte Carlo
#' permutation tests where permutations are sampled with replacement,
#' as derived by Phipson & Smyth (2010). Using the uncorrected estimator
#' `k / n` produces p-values that are systematically understated by
#' approximately `1 / n`, which becomes serious in multiple testing
#' contexts and can produce impossible p-values of exactly zero.
#'
#' Non-finite values in `null_dist` are dropped before computation.
#'
#' The minimum achievable p-value is `1 / (length(null_dist) + 1)`, so
#' Monte Carlo p-values are always bounded strictly above zero.
#'
#' @param observed Numeric scalar. Observed test statistic, e.g. the
#'   F-statistic returned by [test_decay_slopes()] at an anchor of
#'   interest.
#' @param null_dist Numeric vector. Null distribution of the same test
#'   statistic, e.g. as returned by [null_distribution_F()]. Non-finite
#'   entries are dropped.
#'
#' @return A numeric scalar between `1 / (n + 1)` and 1, where `n` is
#'   the number of finite values in `null_dist`.
#'
#' @references
#' Phipson B, Smyth GK (2010). Permutation p-values should never be
#' zero: calculating exact p-values when permutations are randomly
#' drawn. *Statistical Applications in Genetics and Molecular Biology*
#' 9(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' @seealso [null_distribution_F()], [test_decay_significance()].
#'
#' @examples
#' set.seed(1)
#' null_F <- rf(500, df1 = 3, df2 = 50)
#' empirical_pvalue(observed = 5, null_dist = null_F)
#' empirical_pvalue(observed = 100, null_dist = null_F)  # hits the floor
#'
#' @export
empirical_pvalue <- function(observed, null_dist) {
    if (!is.numeric(observed) || length(observed) != 1L ||
        !is.finite(observed)) {
        stop("`observed` must be a single finite numeric value.")
    }
    null_dist <- null_dist[is.finite(null_dist)]
    if (length(null_dist) == 0L) {
        stop("`null_dist` contains no finite values.")
    }
    (sum(null_dist >= observed) + 1L) / (length(null_dist) + 1L)
}


#' Test whether a specific anchor position shows an unusual group effect
#'
#' End-to-end wrapper for anchor-specificity testing. Fits the LMM at
#' the anchor position of interest to obtain an observed F-statistic
#' for the distance-by-group interaction, then builds a null
#' distribution of F-statistics at `n_random` random anchor positions
#' and returns an empirical p-value.
#'
#' @details
#'
#' This is the primary tool for answering the scientific question:
#' "does coverage decay around this specific anchor differ between
#' treatment groups by more than would be expected at a random
#' genomic position?"
#'
#' The empirical p-value complements, rather than replaces, the LMM
#' p-value at the anchor. The two answer different questions:
#'
#' - The LMM p-value asks: "given the observed data, how surprising
#'   is the observed group-by-distance interaction under the null
#'   that groups share a common slope?"
#' - The empirical p-value asks: "how surprising is the observed
#'   F-statistic under the null that there is nothing special about
#'   this anchor position, i.e. we would see a similar effect at any
#'   random genomic position?"
#'
#' A small LMM p-value combined with a large empirical p-value
#' indicates a real group effect that is not specific to the anchor.
#' A small empirical p-value with a small LMM p-value is the
#' publication-grade result: the effect is real *and* localized to
#' the anchor.
#'
#' Runtime is dominated by the `n_random` LMM refits inside the null
#' distribution, so use `BPPARAM = BiocParallel::MulticoreParam()` for
#' meaningful `n_random`.
#'
#' @param combined A data.frame or tibble from [load_all_samples()].
#' @param anchor_position Non-negative numeric scalar. Anchor
#'   coordinate of scientific interest, in base pairs (e.g. OriLyt
#'   coordinate).
#' @param genome_length Positive integer. Reference genome length,
#'   used to bound random anchor sampling.
#' @param n_random Positive integer. Number of random anchor positions
#'   used to build the null. Default 500.
#' @param window Optional positive numeric scalar. Restrict the fit to
#'   bins within `window` bp of the anchor. Passed through to both the
#'   observed fit and every random-anchor fit so the two are
#'   comparable.
#' @param exclude_radius Non-negative numeric scalar. If positive,
#'   exclude a band of width `2 * exclude_radius` centered on
#'   `anchor_position` from the random-anchor pool. Default 0 (no
#'   exclusion). Set to `window` (or larger) to guarantee no overlap
#'   between random-anchor windows and the observed-anchor window.
#' @param seed Optional integer. Random seed for reproducibility of
#'   the sampled random anchor positions.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` instance
#'   controlling parallelism. Default [BiocParallel::SerialParam()].
#'
#' @return An S3 object of class `"decay_significance"`, a named list
#'   with elements:
#'   \describe{
#'     \item{`anchor_test`}{The `decay_test` object at
#'       `anchor_position`, as returned by [test_decay_slopes()].}
#'     \item{`null_F`}{Numeric vector of F-statistics from the random
#'       anchor positions (successful fits only).}
#'     \item{`empirical_p`}{Numeric scalar. Empirical p-value for
#'       the observed F under the null distribution.}
#'     \item{`anchor_position`}{Anchor coordinate echoed for
#'       provenance.}
#'     \item{`genome_length`, `n_random`, `window`}{Inputs echoed for
#'       provenance.}
#'   }
#'
#' @references
#' Phipson B, Smyth GK (2010). Permutation p-values should never be
#' zero: calculating exact p-values when permutations are randomly
#' drawn. *Statistical Applications in Genetics and Molecular Biology*
#' 9(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' @seealso [test_decay_slopes()], [null_distribution_F()],
#'   [empirical_pvalue()].
#'
#' @examples
#' \donttest{
#' ## Illustration only; this call is expensive because it refits the
#' ## LMM at n_random + 1 anchor positions. For real analyses, use
#' ## n_random = 500 or higher and pass a parallel BPPARAM.
#' # combined <- load_all_samples(meta, base_dir, bin_size = 1000,
#' #                              genome_length = 238080)
#' # res <- test_decay_significance(
#' #     combined,
#' #     anchor_position = 93000,       # HCMV OriLyt (verify vs. your reference)
#' #     genome_length   = 238080,
#' #     n_random        = 500,
#' #     window          = 50000,
#' #     exclude_radius  = 50000,
#' #     seed            = 2026,
#' #     BPPARAM         = BiocParallel::MulticoreParam(workers = 8)
#' # )
#' # res$empirical_p
#' }
#'
#' @importFrom BiocParallel SerialParam
#' @export
test_decay_significance <- function(combined,
                                    anchor_position,
                                    genome_length,
                                    n_random = 500L,
                                    window = NULL,
                                    exclude_radius = 0,
                                    seed = NULL,
                                    BPPARAM = BiocParallel::SerialParam()) {
    .validate_position(anchor_position, "anchor_position",
        genome_length = genome_length)
    if (!is.numeric(exclude_radius) || length(exclude_radius) != 1L ||
        !is.finite(exclude_radius) || exclude_radius < 0) {
        stop("`exclude_radius` must be a single non-negative number.")
    }

    anchor_test <- test_decay_slopes(
        combined,
        anchor_position = anchor_position,
        window          = window,
        quiet           = TRUE
    )

    exclude <- if (exclude_radius > 0) {
        c(max(0, anchor_position - exclude_radius),
            min(genome_length, anchor_position + exclude_radius))
    } else NULL

    null_F <- null_distribution_F(
        combined,
        genome_length = genome_length,
        n_random      = n_random,
        window        = window,
        exclude       = exclude,
        seed          = seed,
        BPPARAM       = BPPARAM
    )

    p_emp <- empirical_pvalue(anchor_test$F_stat, null_F)

    structure(
        list(
            anchor_test     = anchor_test,
            null_F          = null_F,
            empirical_p     = p_emp,
            anchor_position = anchor_position,
            genome_length   = genome_length,
            n_random        = n_random,
            window          = window
        ),
        class = "decay_significance"
    )
}


#' Print method for decay_significance objects
#'
#' Prints a compact summary of a [test_decay_significance()] result,
#' showing the observed LMM F-statistic and p-value at the anchor,
#' the number of successful null fits, and the empirical p-value.
#'
#' @param x A `decay_significance` object.
#' @param digits Positive integer. Significant digits for numeric
#'   formatting. Default 4.
#' @param ... Ignored.
#'
#' @return `x` invisibly.
#'
#' @export
print.decay_significance <- function(x, digits = 4L, ...) {
    cat("Coverage decay significance test\n")
    cat("--------------------------------\n")
    cat(sprintf("Anchor position   : %s\n", format(x$anchor_position)))
    cat(sprintf("Window            : %s\n",
        if (is.null(x$window)) "full genome"
        else sprintf("+/- %g bp", x$window)))
    cat(sprintf("Observed F        : %s\n",
        format(round(x$anchor_test$F_stat, digits))))
    cat(sprintf("Anchor LMM p      : %s\n",
        format.pval(x$anchor_test$p_value, digits = digits)))
    cat(sprintf("Null F fits       : %d of %d requested\n",
        length(x$null_F), x$n_random))
    cat(sprintf("Empirical p-value : %s\n",
        format.pval(x$empirical_p, digits = digits)))
    invisible(x)
}
