## =========================================================================
## Sensitivity modifiers S1 to S5 (PLAN.md 2.3). SCAFFOLDING ONLY.
##
## A2-1 defines the hook signatures that `generate_replicate()` calls so the
## baseline generator has clean extension points. Full implementations are
## A2-3. A hook runs only when the scenario's `sens` list names it, so the
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
#' @param e Numeric matrix, atoms x samples, independent N(0, sigma_atom^2).
#' @param sc A `cda_scenario`; `sc$sens$S1a$phi` is the AR(1) coefficient.
#' @return Matrix of the same shape with AR(1) columns.
sens_s1a_noise <- function(e, sc) {
    .not_implemented("S1a")
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
    .not_implemented("S1b")
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
    .not_implemented("S2")
}

#' S3. Heteroscedasticity across groups.
#'
#' Group-specific random-slope SD (`tau_g = cv * s_g`) and/or a residual SD
#' ratio between groups.
#'
#' @param samples Tibble of per-sample truth with `group` and `s_ref_g`.
#' @param sc A `cda_scenario`; `sc$sens$S3` holds `cv` and `sigma_ratio`.
#' @return `samples` with columns `tau_i` and `sigma_i` (per-sample SDs).
sens_s3_scales <- function(samples, sc) {
    .not_implemented("S3")
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
