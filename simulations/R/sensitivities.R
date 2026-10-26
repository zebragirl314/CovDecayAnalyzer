## =========================================================================
## Sensitivity modifiers S1 to S5 (PLAN.md 2.3).
##
## A2-1 defined the hook signatures that `generate_replicate()` calls so the
## baseline generator has clean extension points. S1a, S1b, S2 and S3 are
## implemented in A2-4; S4 (low depth) and S5 (anchor misplacement) remain
## stubs (A2-4's "anchor at extremes" needs no hook: it is `anchor_frac` 0 or
## 1). A hook runs only when the scenario's `sens` list names it, so the
## baseline path never reaches these stubs. Each sensitivity is a one-factor
## modification of the baseline and consumes random numbers only *after*
## the baseline draws, so baseline and sensitivity datasets built from the
## same seed share their baseline components (paired deltas).
##
## `sc$sens` is a named list, for example
##   list(S1a = list(phi = 0.6), S2 = list(frac = 0.02, run_kb = 5, ...)).
## S6 (end anchors) needs no hook: it is `anchor_frac` = 0 or 1.
## =========================================================================

sens_active <- function(sc, name) {
    !is.null(sc$sens[[name]])
}

.not_implemented <- function(name) {
    stop(sprintf("Sensitivity %s is scaffolded only; implemented in A2-3.",
        name), call. = FALSE)
}

#' S1a. Stationary residual autocorrelation.
#'
#' Replaces independent atom noise with AR(1) noise along the genome,
#' coefficient `phi`, marginal SD held at `sigma_atom`.
#'
#' The independent draws are reused as the innovations, so the AR(1) dataset
#' shares its baseline draws with the baseline scenario: `x_1 = e_1`,
#' `x_t = phi * x_{t-1} + sqrt(1 - phi^2) * e_t`.
#'
#' @param e Numeric matrix, atoms x samples, independent N(0, sigma_i^2).
#' @param sc A `cda_scenario`; `sc$sens$S1a$phi` is the AR(1) coefficient.
#' @return Matrix of the same shape with AR(1) columns.
sens_s1a_noise <- function(e, sc) {
    phi <- sc$sens$S1a$phi
    if (!is.numeric(phi) || length(phi) != 1L || abs(phi) >= 1) {
        stop("`sens$S1a$phi` must be a single number in (-1, 1).")
    }
    innov <- sqrt(1 - phi^2)
    x <- e
    for (t in seq_len(nrow(e))[-1L]) {
        x[t, ] <- phi * x[t - 1L, ] + innov * e[t, ]
    }
    x
}

#' S1b. Systematic misfit of the mean curve.
#'
#' Replaces the log-linear mean with `f + (1 - f) * exp(-s_i * d)` (floor
#' `f`, optionally sample-specific). Under S1b `Delta` is a rate difference
#' in the true curve; the LMM slope is a biased summary of it.
#'
#' @param d Numeric vector of atom distances from the true anchor.
#' @param s_i Numeric vector of per-sample decay rates.
#' @param sc A `cda_scenario`; `sc$sens$S1b$floor` is `f`.
#' @return Matrix (atoms x samples) of log mean coverage, before depth and
#'   noise.
sens_s1b_curve <- function(d, s_i, sc) {
    f <- sc$sens$S1b$floor
    if (!is.numeric(f) || length(f) != 1L || f < 0 || f >= 1) {
        stop("`sens$S1b$floor` must be a single number in [0, 1).")
    }
    log(f + (1 - f) * exp(-outer(d, s_i)))
}

#' S2. Zero-coverage dropout runs.
#'
#' Sets runs of atoms to exactly zero (per-sample Poisson run starts,
#' geometric run lengths) before binning.
#'
#' @param cov Numeric matrix, atoms x samples, mean coverage.
#' @param sc A `cda_scenario`; `sc$sens$S2` holds fraction, mean run length,
#'   placement, and affected samples.
#' @return Matrix of the same shape with zeroed runs.
sens_s2_dropout <- function(cov, sc) {
    frac <- sc$sens$S2$frac
    run_kb <- sc$sens$S2$run_kb
    if (!is.numeric(frac) || frac < 0 || frac >= 1 ||
        !is.numeric(run_kb) || run_kb < 1) {
        stop("`sens$S2` needs `frac` in [0, 1) and `run_kb` >= 1.")
    }
    n_atoms <- nrow(cov)
    mean_len <- run_kb * 1000 / sc$atom_size  # mean run length in atoms
    ## Starts per atom chosen so the expected zeroed fraction is `frac`
    ## (overlaps and truncation at the genome ends make it slightly lower).
    lambda <- frac / mean_len
    for (i in seq_len(ncol(cov))) {
        k <- stats::rpois(1L, lambda * n_atoms)
        if (k == 0L) next
        starts <- sample.int(n_atoms, k, replace = TRUE)
        lens <- stats::rgeom(k, 1 / mean_len) + 1L
        for (j in seq_len(k)) {
            cov[starts[j]:min(n_atoms, starts[j] + lens[j] - 1L), i] <- 0
        }
    }
    cov
}

#' S3. Heteroscedasticity across groups.
#'
#' Group-specific random-slope CV (`cv_g`) and/or a residual SD
#' ratio between groups.
#'
#' The last group (the one carrying the effect) has its relative slope SD
#' multiplied by `cv_ratio` and its atom residual SD by `sigma_ratio`; the
#' other groups keep the scenario's `cv` and `sigma_atom`. Both default to 1.
#' The random draws are scaled, not replaced, so the dataset pairs with the
#' baseline scenario.
#'
#' @param samples Tibble of per-sample truth with `group`, `cv_i`, `sigma_i`.
#' @param sc A `cda_scenario`; `sc$sens$S3` holds `cv_ratio` and
#'   `sigma_ratio`.
#' @return `samples` with columns `cv_i` (relative slope SD) and `sigma_i`.
sens_s3_scales <- function(samples, sc) {
    cv_ratio <- sc$sens$S3$cv_ratio %||% 1
    sigma_ratio <- sc$sens$S3$sigma_ratio %||% 1
    if (!is.numeric(cv_ratio) || !is.numeric(sigma_ratio) ||
        cv_ratio <= 0 || sigma_ratio <= 0) {
        stop("`sens$S3` ratios must be positive numbers.")
    }
    last <- samples$group == levels(samples$group)[nlevels(samples$group)]
    samples$cv_i[last] <- samples$cv_i[last] * cv_ratio
    samples$sigma_i[last] <- samples$sigma_i[last] * sigma_ratio
    samples
}

#' S4. Low depth with Poisson counting noise.
#'
#' Replaces expected atom coverage with a Poisson count layer (the median
#' depth itself is set through `sc$depth_median`).
#'
#' @param cov Numeric matrix, atoms x samples, expected mean coverage.
#' @param sc A `cda_scenario`.
#' @return Matrix of the same shape with counting noise.
sens_s4_counts <- function(cov, sc) {
    .not_implemented("S4")
}

#' S5. Anchor misplacement.
#'
#' The analysis anchor is offset from the true anchor by
#' `sc$sens$S5$offset` bp.
#'
#' @param sc A `cda_scenario`.
#' @return Numeric scalar, the anchor used for analysis.
sens_s5_analysis_anchor <- function(sc) {
    .not_implemented("S5")
}
