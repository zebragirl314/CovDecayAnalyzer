## =========================================================================
## Linear mixed-effects test of group differences in coverage decay.
##
## The model:
##   log_norm_cov ~ distance_from_anchor * group +
##                  (1 + distance_from_anchor | sample)
##
## The null hypothesis of interest is that the distance-by-group interaction
## coefficient(s) are zero, i.e. all groups share a common decay slope. This
## is tested via a Kenward-Roger F test, which provides small-sample-
## corrected denominator degrees of freedom appropriate for the low
## per-group replication (n = 2 to 5) typical of these experiments.
##
## When the interaction is significant, pairwise group contrasts on the
## decay slope are obtained from emmeans.
##
## Numerical stability note: raw distances in base pairs can be very large
## (100 kb or more for viral genomes) while per-bp slopes are correspondingly
## tiny (~1e-5). The resulting design matrix condition number is far beyond
## double precision and Kenward-Roger's variance adjustment step fails with
## "system is computationally singular". We rescale distance internally to
## put its maximum in the 1-to-10 range and rescale reported slope estimates
## back to per-bp on output, so the user-facing API is unchanged.
## =========================================================================


#' Test group differences in coverage decay using a linear mixed model
#'
#' Fits a linear mixed-effects model to all bins across all samples
#' simultaneously and tests whether the rate of coverage decay with
#' distance from the anchor differs between treatment groups. This is
#' the recommended formal test in the CovDecayAnalyzer workflow.
#'
#' @section Statistical model:
#'
#' The fitted model is
#'
#' \deqn{y_{ij} = (\beta_0 + b_{0i}) + (\beta_1 + b_{1i}) \cdot d_{ij} +
#'                \sum_{g=2}^{G} \gamma_g \cdot I(group_i = g) +
#'                \sum_{g=2}^{G} \delta_g \cdot d_{ij} \cdot I(group_i = g) +
#'                \varepsilon_{ij}}{
#' y_ij = (b0 + u0_i) + (b1 + u1_i) * d_ij + sum(gamma_g * I(group_i = g))
#'        + sum(delta_g * d_ij * I(group_i = g)) + e_ij}
#'
#' where `y_ij` is `log_norm_cov` for bin `j` of sample `i`, `d_ij` is
#' `distance_from_anchor`, `beta_0` and `beta_1` are the fixed intercept
#' and slope for the reference group, `u_0i` and `u_1i` are per-sample
#' random intercept and slope deviations, and `gamma_g`, `delta_g` are
#' fixed group main-effects and group-by-distance interaction terms.
#'
#' The random-effects structure `(1 + distance_from_anchor | sample)`
#' allows each sample to have its own intercept and its own slope. This
#' properly represents the fact that different samples have different
#' overall coverage levels and possibly different local decay rates
#' even within the same group.
#'
#' @section Test and small-sample correction:
#'
#' The scientific hypothesis of "do groups differ in decay rate?"
#' corresponds to the joint hypothesis
#' `H0: delta_2 = delta_3 = ... = delta_G = 0`. This is tested as a
#' Type III F test on the `distance_from_anchor:group` interaction
#' using the Kenward-Roger denominator-degree-of-freedom approximation
#' (Kenward & Roger 1997), as implemented in `lmerTest` (Kuznetsova,
#' Brockhoff & Christensen 2017) on top of `lme4` (Bates et al. 2015).
#'
#' The Kenward-Roger correction is specifically designed for the
#' low-replication regime typical of these experiments (n = 2 to 5
#' samples per group). It adjusts both the fixed-effect variance
#' estimate and the denominator degrees of freedom to control Type I
#' error at the nominal level. Naive Wald or Satterthwaite tests can
#' be liberal in this regime; Kenward-Roger is the standard
#' recommendation.
#'
#' @section Pairwise contrasts:
#'
#' When the overall interaction is significant, per-group decay slopes
#' and Tukey-adjusted pairwise contrasts are extracted using
#' [emmeans::emtrends()] (Lenth 2024). This returns each group's
#' estimated slope with a confidence interval, and each pairwise slope
#' difference with a p-value adjusted for the number of contrasts.
#'
#' @section Numerical scaling:
#'
#' Distance from the anchor is internally rescaled to put its maximum
#' in the 1-to-10 range before fitting the LMM. This keeps the design
#' matrix well-conditioned and lets the Kenward-Roger variance
#' adjustment converge on inputs spanning tens or hundreds of kilobases,
#' where the raw-bp representation is numerically singular. Slope
#' estimates and their standard errors are converted back to per-bp
#' units on the way out, so the returned `group_slopes` and `pairwise`
#' tables report slopes in the same units the user expects
#' (log_norm_cov per base pair). The F-statistic and p-value are
#' invariant to this rescaling.
#'
#' @param combined A data.frame or tibble from [load_all_samples()].
#'   Must contain columns `sample`, `group`, `bin_mid`, and
#'   `log_norm_cov`. If a `distance_from_anchor` column is already
#'   present, `anchor_position` is ignored.
#' @param anchor_position Non-negative numeric scalar. Reference
#'   coordinate for distance calculation. Ignored if `combined`
#'   already contains a `distance_from_anchor` column.
#' @param window Optional positive numeric scalar. If supplied, only
#'   bins within `window` base pairs of the anchor are used in the
#'   fit.
#' @param reml Logical. If `TRUE` (default), fit by restricted maximum
#'   likelihood (REML). REML is required for correct Kenward-Roger
#'   inference and should not be changed unless you have a specific
#'   reason.
#' @param quiet Logical. If `TRUE` (default), suppress convergence and
#'   singular-fit messages from `lme4`. Fitting errors are still raised
#'   to the caller.
#'
#' @return An S3 object of class `"decay_test"`, a named list with
#'   elements:
#'   \describe{
#'     \item{`model`}{Fitted `lmerModLmerTest` object. Note that the
#'       fit's `distance_scaled` variable is internally rescaled for
#'       numerical stability; see Numerical scaling.}
#'     \item{`anova_table`}{Kenward-Roger ANOVA table for the fixed
#'       effects.}
#'     \item{`F_stat`}{Numeric scalar. F-statistic for the
#'       distance-by-group interaction.}
#'     \item{`df1`, `df2`}{Numerator and Kenward-Roger denominator
#'       degrees of freedom for the interaction test.}
#'     \item{`p_value`}{Numeric scalar. P-value for the interaction
#'       test.}
#'     \item{`group_slopes`}{Tibble of estimated per-group slopes
#'       (per base pair) with confidence intervals.}
#'     \item{`pairwise`}{Tibble of pairwise group slope contrasts
#'       (per base pair), Tukey-adjusted.}
#'     \item{`n_samples`, `n_bins`}{Fit summary statistics.}
#'     \item{`anchor_position`, `window`}{Inputs echoed for
#'       provenance.}
#'     \item{`distance_scale`}{Internal scaling factor used to
#'       rescale distance before fitting.}
#'   }
#'
#' @references
#' Bates D, Machler M, Bolker B, Walker S (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software*
#' 67(1):1-48. \doi{10.18637/jss.v067.i01}
#'
#' Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package:
#' Tests in linear mixed effects models. *Journal of Statistical Software*
#' 82(13):1-26. \doi{10.18637/jss.v082.i13}
#'
#' Kenward MG, Roger JH (1997). Small sample inference for fixed effects
#' from restricted maximum likelihood. *Biometrics* 53(3):983-997.
#' \doi{10.2307/2533558}
#'
#' Lenth RV (2024). emmeans: Estimated Marginal Means, aka Least-Squares
#' Means. R package.
#' \url{https://CRAN.R-project.org/package=emmeans}
#'
#' @seealso
#' [fit_decay_slopes()] for a per-sample diagnostic view;
#' [test_decay_significance()] for the wrapper that combines this test
#' with a random-anchor null distribution.
#'
#' @examples
#' \donttest{
#' extdata_dir <- system.file("extdata", package = "CovDecayAnalyzer")
#' meta <- read.csv(
#'     file.path(extdata_dir, "sample_metadata.csv"),
#'     stringsAsFactors = FALSE
#' )
#' meta$coverage_file <- breseq_coverage_paths(meta$sample, extdata_dir)
#' combined <- load_all_samples(
#'     sample_metadata = meta,
#'     bin_size        = 100,
#'     genome_length   = 2000
#' )
#' res <- test_decay_slopes(combined, anchor_position = 1000)
#' res
#' res$pairwise
#' }
#'
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
#' @importFrom stats anova
#' @importFrom emmeans emtrends
#' @importFrom graphics pairs
#' @importFrom tibble as_tibble
#' @export
test_decay_slopes <- function(combined,
                            anchor_position = NULL,
                            window = NULL,
                            reml = TRUE,
                            quiet = TRUE) {
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
    if (!is.logical(reml) || length(reml) != 1L) {
        stop("`reml` must be TRUE or FALSE.")
    }

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

    combined <- combined[is.finite(combined$log_norm_cov) &
        is.finite(combined$distance_from_anchor), , drop = FALSE]
    combined$group <- factor(combined$group)

    n_groups <- length(levels(combined$group))
    if (n_groups < 2L) {
        stop("Need at least two groups to test group differences.")
    }
    n_samples <- length(unique(combined$sample))
    if (n_samples < 3L) {
        stop("Need at least three samples for random-effects fitting.")
    }

    ## ---- Numerical rescaling of distance for LMM stability ------------
    ## Choose a scale so max(distance) lands in the 1-to-10 range. This
    ## keeps the design matrix condition number in a range where
    ## Kenward-Roger's variance adjustment can converge without
    ## reporting "system is computationally singular".
    max_dist <- max(combined$distance_from_anchor, na.rm = TRUE)
    dist_scale <- if (max_dist > 10) max_dist / 10 else 1
    combined$distance_scaled <- combined$distance_from_anchor / dist_scale

    .maybe_quiet <- function(expr) {
        if (isTRUE(quiet)) {
            suppressMessages(suppressWarnings(expr))
        } else {
            expr
        }
    }

    fit <- .maybe_quiet(
        lmerTest::lmer(
            log_norm_cov ~ distance_scaled * group +
                (1 + distance_scaled | sample),
            data = combined,
            REML = reml,
            control = lme4::lmerControl(
                check.conv.singular = "ignore",
                check.conv.grad     = "ignore",
                check.conv.hess     = "ignore"
            )
        )
    )

    ## Kenward-Roger ANOVA. F, df, and p are scale-invariant so no
    ## adjustment needed on this table.
    anova_tbl <- .maybe_quiet(
        stats::anova(fit, ddf = "Kenward-Roger")
    )
    interaction_row <- grep("distance_scaled:group",
        rownames(anova_tbl), value = TRUE)
    if (length(interaction_row) != 1L) {
        stop("Could not locate `distance_scaled:group` row in the ",
            "Kenward-Roger ANOVA table.")
    }
    F_stat <- anova_tbl[interaction_row, "F value"]
    df1    <- anova_tbl[interaction_row, "NumDF"]
    df2    <- anova_tbl[interaction_row, "DenDF"]
    p_val  <- anova_tbl[interaction_row, "Pr(>F)"]

    ## Rename ANOVA rownames to user-facing variable name.
    rownames(anova_tbl) <- gsub("distance_scaled",
        "distance_from_anchor", rownames(anova_tbl), fixed = TRUE)

    ## Per-group slope estimates and pairwise contrasts from emmeans.
    ## Note: `distance_scaled` is the variable emmeans knows about.
    ## The slope estimates come back in units of log_norm_cov per
    ## scaled-distance-unit. To report per bp, divide by dist_scale.
    emt <- .maybe_quiet(
        emmeans::emtrends(fit, ~ group, var = "distance_scaled")
    )
    group_slopes <- as.data.frame(emt)
    group_slopes <- .rescale_slope_frame(group_slopes, dist_scale,
        old_prefix = "distance_scaled",
        new_prefix = "distance_from_anchor")

    pairwise <- as.data.frame(
        .maybe_quiet(pairs(emt, adjust = "tukey"))
    )
    pairwise <- .rescale_slope_frame(pairwise, dist_scale)

    structure(
        list(
            model           = fit,
            anova_table     = anova_tbl,
            F_stat          = as.numeric(F_stat),
            df1             = as.numeric(df1),
            df2             = as.numeric(df2),
            p_value         = as.numeric(p_val),
            group_slopes    = tibble::as_tibble(group_slopes),
            pairwise        = tibble::as_tibble(pairwise),
            n_samples       = n_samples,
            n_bins          = nrow(combined),
            anchor_position = anchor_position,
            window          = window,
            distance_scale  = dist_scale
        ),
        class = "decay_test"
    )
}


## Internal: rescale slope-like columns of an emmeans summary
## data.frame back to per-bp units. Optionally renames a column prefix
## (used by emtrends output where the trend column is named
## "distance_scaled.trend").
##
## Columns rescaled: any *.trend column, plus estimate, SE, lower.CL,
## upper.CL, asymp.LCL, asymp.UCL. The t.ratio, df, and p.value columns
## are scale-invariant and left as-is.
#' @keywords internal
#' @noRd
.rescale_slope_frame <- function(df, dist_scale,
                                 old_prefix = NULL,
                                 new_prefix = NULL) {
    slope_cols <- c(
        grep("\\.trend$", names(df), value = TRUE),
        intersect(names(df),
            c("estimate", "SE", "lower.CL", "upper.CL",
                "asymp.LCL", "asymp.UCL"))
    )
    for (cc in slope_cols) {
        df[[cc]] <- df[[cc]] / dist_scale
    }
    if (!is.null(old_prefix) && !is.null(new_prefix)) {
        names(df) <- gsub(old_prefix, new_prefix, names(df), fixed = TRUE)
    }
    df
}


#' Print method for decay_test objects
#'
#' Prints a compact summary of a [test_decay_slopes()] result.
#'
#' @param x A `decay_test` object.
#' @param digits Positive integer. Significant digits for numeric
#'   formatting. Default 4.
#' @param ... Ignored.
#'
#' @return `x` invisibly.
#'
#' @export
print.decay_test <- function(x, digits = 4L, ...) {
    cat("Coverage decay LMM test\n")
    cat("-----------------------\n")
    cat(sprintf("Anchor position : %s\n",
        if (is.null(x$anchor_position)) "<supplied via combined>"
        else format(x$anchor_position)))
    cat(sprintf("Window          : %s\n",
        if (is.null(x$window)) "full genome"
        else sprintf("+/- %g bp", x$window)))
    cat(sprintf("Samples fit     : %d\n", x$n_samples))
    cat(sprintf("Bins fit        : %d\n", x$n_bins))
    cat("\nFixed-effect interaction (distance x group)\n")
    cat(sprintf("  F(%g, %.2f) = %s   p = %s\n",
        x$df1, x$df2,
        format(round(x$F_stat, digits)),
        format.pval(x$p_value, digits = digits)))
    cat("\nPer-group decay slopes (emtrends, per bp)\n")
    print(x$group_slopes, n = Inf)
    cat("\nPairwise group contrasts (Tukey-adjusted, per bp)\n")
    print(x$pairwise, n = Inf)
    invisible(x)
}
