## =========================================================================
## Visualization functions for coverage decay analysis.
##
## All functions return a ggplot object rather than printing it, so the
## caller can modify layers, save with ggsave, or combine with other
## plots via patchwork or cowplot.
## =========================================================================


## Okabe-Ito colorblind-safe palette (Okabe & Ito 2008), the default for
## categorical colors in every plot function.
.okabe_ito <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
)

## Palette names accepted by ggplot2's Brewer scales. Names in this set
## keep routing through `scale_*_brewer()` for back-compatibility.
.brewer_palette_names <- c(
    "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd",
    "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu",
    "YlOrBr", "YlOrRd",
    "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn",
    "Spectral",
    "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2",
    "Set3"
)

## Discrete color or fill scale for `n_levels` groups. A Brewer palette
## name uses the Brewer scale; anything else uses the Okabe-Ito vector.
## `drop = FALSE` keeps each group on its own palette position when a
## plot shows only some of the levels (see `plot_sample()`).
#' @importFrom ggplot2 scale_color_brewer scale_fill_brewer scale_color_manual scale_fill_manual
#' @noRd
.palette_scale <- function(palette, n_levels, aesthetic = c("color", "fill")) {
    aesthetic <- match.arg(aesthetic)
    if (!is.character(palette) || length(palette) != 1L || is.na(palette)) {
        stop("`palette` must be a single character string.")
    }
    if (palette %in% .brewer_palette_names) {
        return(if (aesthetic == "color") {
            ggplot2::scale_color_brewer(palette = palette, drop = FALSE)
        } else {
            ggplot2::scale_fill_brewer(palette = palette, drop = FALSE)
        })
    }
    if (palette != "okabe_ito") {
        warning(sprintf(
            "Unknown palette '%s'; using the Okabe-Ito palette.", palette
        ))
    }
    if (n_levels > length(.okabe_ito)) {
        stop(sprintf(
            "The Okabe-Ito palette has %d colors but there are %d groups; ",
            length(.okabe_ito), n_levels
        ), "supply a Brewer palette name via `palette`.")
    }
    values <- .okabe_ito[seq_len(n_levels)]
    if (aesthetic == "color") {
        ggplot2::scale_color_manual(values = values, drop = FALSE)
    } else {
        ggplot2::scale_fill_manual(values = values, drop = FALSE)
    }
}


#' Plot normalized coverage vs distance from anchor, colored by group
#'
#' Produces a ggplot showing each sample's per-bin normalized coverage
#' as a function of distance from the anchor position, with a
#' loess-smoothed trend line per treatment group. This is the primary
#' visual diagnostic for coverage decay analyses.
#'
#' @details
#'
#' The plot layers, in order:
#'
#' - Optional per-bin points (`geom_point`) with low alpha, so
#'   individual observations are visible but the group trends dominate.
#' - A loess-smoothed line per group (`geom_smooth`) with 95% confidence
#'   ribbon. Loess is used rather than the log-linear fit from the
#'   formal test because visualization benefits from a nonparametric
#'   smoother that can reveal deviations from linearity.
#'
#' The loess `span` parameter controls the degree of smoothing: smaller
#' values follow the data more closely and can reveal local structure,
#' larger values produce smoother curves. The default of 0.3 works well
#' for a full-genome view at 1 kb binning; use 0.5 or higher for very
#' small windows where 0.3 would be unstable.
#'
#' If `combined` already contains a `distance_from_anchor` column, it is
#' used as-is. Otherwise, `anchor_position` is required and distances
#' are computed via [add_distance_from_anchor()].
#'
#' @param combined A data.frame or tibble from [load_all_samples()].
#'   Must contain columns `sample`, `group`, `bin_mid`, and `norm_cov`.
#'   May optionally contain `distance_from_anchor`.
#' @param anchor_position Non-negative numeric scalar. Reference
#'   coordinate for distance calculation. Ignored if `combined`
#'   already contains `distance_from_anchor`.
#' @param window Optional positive numeric scalar. If supplied,
#'   restricts the plot's x-axis (via `coord_cartesian`) to
#'   `[0, window]`. Data outside the window are still fit by the loess
#'   smoother but are not shown; this is preferable to filtering the
#'   data before smoothing, which would truncate the smoother's
#'   support and introduce edge artifacts.
#' @param show_points Logical. If `TRUE` (default), overlay individual
#'   per-bin points on top of the group trends.
#' @param loess_span Numeric in `(0, 1]`. Passed to `geom_smooth(span =
#'   loess_span, method = "loess")`. Default 0.3.
#' @param palette Character scalar. `"okabe_ito"` (default) uses the
#'   Okabe-Ito colorblind-safe palette, which has eight colors, so at
#'   most eight groups. A Brewer palette name (e.g. `"Set1"`) uses
#'   [ggplot2::scale_color_brewer()] instead. Any other string falls back
#'   to Okabe-Ito with a warning.
#'
#' @return A `ggplot` object. Not printed automatically; call `print()`
#'   or assign into a knitr code chunk, or save with
#'   [ggplot2::ggsave()].
#'
#' @seealso [plot_decay_slopes()], [plot_null_distribution()].
#'
#' @examples
#' ## End-to-end example on the shipped simulated data. The WT+GCV
#' ## samples decay six times faster than the WT-GCV samples, so the
#' ## two loess curves should be clearly separated.
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
#' plot_coverage_decay(combined, anchor_position = 1000, window = 1000)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth coord_cartesian labs theme_minimal theme element_text
#' @export
plot_coverage_decay <- function(combined,
                                anchor_position = NULL,
                                window = NULL,
                                show_points = TRUE,
                                loess_span = 0.3,
                                palette = "okabe_ito") {
    if (!is.data.frame(combined)) {
        stop("`combined` must be a data.frame or tibble.")
    }
    required <- c("sample", "group", "bin_mid", "norm_cov")
    missing_cols <- setdiff(required, names(combined))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "`combined` is missing required column(s): %s",
            paste(missing_cols, collapse = ", ")
        ))
    }
    if (!is.numeric(loess_span) || length(loess_span) != 1L ||
        loess_span <= 0 || loess_span > 1) {
        stop("`loess_span` must be a single number in (0, 1].")
    }
    if (!"distance_from_anchor" %in% names(combined)) {
        if (is.null(anchor_position)) {
            stop("Provide `anchor_position` or a `distance_from_anchor` ",
                "column in `combined`.")
        }
        combined <- add_distance_from_anchor(combined, anchor_position)
    }

    combined$group <- factor(combined$group)
    group_scale <- .palette_scale(palette, nlevels(combined$group))

    p <- ggplot2::ggplot(
        combined,
        ggplot2::aes(x = distance_from_anchor, y = norm_cov,
            color = group)
    )
    if (isTRUE(show_points)) {
        p <- p + ggplot2::geom_point(alpha = 0.15, size = 0.5)
    }
    p <- p +
        ggplot2::geom_smooth(
            method = "loess", span = loess_span,
            se = TRUE, linewidth = 1.2, formula = y ~ x
        ) +
        group_scale +
        ggplot2::labs(
            title = "Normalized coverage vs distance from anchor",
            x = "Distance from anchor (bp)",
            y = "Normalized coverage (median-scaled)",
            color = "Group"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            plot.title      = ggplot2::element_text(face = "bold"),
            legend.position = "top"
        )
    if (!is.null(window)) {
        if (!is.numeric(window) || length(window) != 1L || window <= 0) {
            stop("`window` must be a single positive number.")
        }
        p <- p + ggplot2::coord_cartesian(xlim = c(0, window))
    }
    p
}


#' Plot per-sample decay slopes by group
#'
#' Draws a boxplot of per-sample decay slopes across treatment groups,
#' overlaid with individual sample points. Takes the tibble output of
#' [fit_decay_slopes()] as input.
#'
#' @details
#'
#' This visualization complements [test_decay_slopes()] by showing the
#' spread and central tendency of per-sample slopes within each group.
#' Consistent slopes within a group and clear separation between groups
#' visually support a formal group difference; overlapping distributions
#' with wide within-group spread suggest the effect is weak or noisy.
#'
#' Note that the boxplot summarizes per-sample slopes as if they were
#' independent point estimates, which is exactly what the formal LMM
#' test *avoids*. Use this plot to understand the data structure and
#' identify outliers, not to draw statistical conclusions.
#'
#' Significance annotations are not added by default. To overlay
#' pairwise comparison bars, extract `x$pairwise` from a `decay_test`
#' object and add a layer such as
#' `ggpubr::stat_pvalue_manual()` manually.
#'
#' @param slopes A tibble as returned by [fit_decay_slopes()]. Must
#'   contain columns `group` and `slope`. `sample` is used for
#'   optional labeling.
#' @param show_jitter Logical. If `TRUE` (default), overlay individual
#'   sample points via `geom_jitter`.
#' @param show_outliers Logical. If `TRUE`, boxplot outliers are drawn
#'   as points. Default `FALSE` (avoids double-plotting when
#'   `show_jitter = TRUE`).
#' @param palette Character scalar. `"okabe_ito"` (default) uses the
#'   Okabe-Ito colorblind-safe palette (at most eight groups). A Brewer
#'   palette name (e.g. `"Set1"`) uses [ggplot2::scale_fill_brewer()]
#'   instead. Any other string falls back to Okabe-Ito with a warning.
#'
#' @return A `ggplot` object.
#'
#' @seealso [fit_decay_slopes()], [plot_coverage_decay()].
#'
#' @examples
#' set.seed(1)
#' fake <- do.call(rbind, lapply(paste0("s", 1:10), function(sname) {
#'     grp <- if (as.integer(sub("s", "", sname)) <= 5) "g1" else "g2"
#'     rate <- if (grp == "g1") -1e-5 else -3e-5
#'     tibble::tibble(
#'         sample  = sname, group = grp,
#'         bin_mid = seq(500, 100000, by = 1000),
#'         log_norm_cov = rate * seq(500, 100000, by = 1000) +
#'             rnorm(100, sd = 0.1)
#'     )
#' }))
#' slopes <- fit_decay_slopes(fake, anchor_position = 0)
#' plot_decay_slopes(slopes)
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme_minimal theme element_text
#' @export
plot_decay_slopes <- function(slopes,
                            show_jitter = TRUE,
                            show_outliers = FALSE,
                            palette = "okabe_ito") {
    if (!is.data.frame(slopes)) {
        stop("`slopes` must be a data.frame or tibble.")
    }
    required <- c("group", "slope")
    missing_cols <- setdiff(required, names(slopes))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "`slopes` is missing required column(s): %s",
            paste(missing_cols, collapse = ", ")
        ))
    }
    slopes$group <- factor(slopes$group)
    group_scale <- .palette_scale(palette, nlevels(slopes$group), "fill")

    outlier_shape <- if (isTRUE(show_outliers)) 19 else NA
    p <- ggplot2::ggplot(
        slopes,
        ggplot2::aes(x = group, y = slope, fill = group)
    ) +
        ggplot2::geom_boxplot(alpha = 0.65, outlier.shape = outlier_shape) +
        group_scale +
        ggplot2::labs(
            title = "Distribution of coverage decay slopes",
            x = "Group",
            y = "Slope (log-linear decay per bp)"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            legend.position = "none",
            plot.title      = ggplot2::element_text(face = "bold")
        )

    if (isTRUE(show_jitter)) {
        p <- p + ggplot2::geom_jitter(width = 0.2, size = 2.5, alpha = 0.8)
    }
    p
}


#' Plot the random-anchor null F distribution with the observed value marked
#'
#' Given a [test_decay_significance()] result, draws a histogram of the
#' null distribution of F-statistics from random anchors and overlays a
#' vertical line at the observed F-statistic at the anchor of interest.
#' The empirical p-value is displayed in the plot title.
#'
#' @details
#'
#' This is the most direct visual representation of the empirical
#' p-value method: the reader can see at a glance whether the observed
#' F falls in the bulk of the null distribution (unremarkable) or in
#' the tail (statistically significant for anchor specificity).
#'
#' @param sig A `decay_significance` object as returned by
#'   [test_decay_significance()].
#' @param bins Positive integer. Number of histogram bins. Default 40.
#' @param observed_color Character scalar. Color of the vertical line
#'   marking the observed F-statistic. Default `"#D55E00"` (Okabe-Ito
#'   vermillion), which stays distinct from the grey histogram.
#'
#' @return A `ggplot` object.
#'
#' @seealso [test_decay_significance()], [null_distribution_F()].
#'
#' @examples
#' set.seed(1)
#' fake_sig <- structure(
#'     list(
#'         anchor_test     = list(F_stat = 8.5, p_value = 0.002),
#'         null_F          = rf(500, df1 = 3, df2 = 40),
#'         empirical_p     = 0.006,
#'         anchor_position = 93000,
#'         genome_length   = 238080,
#'         n_random        = 500,
#'         window          = 50000
#'     ),
#'     class = "decay_significance"
#' )
#' plot_null_distribution(fake_sig)
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_minimal theme element_text
#' @export
plot_null_distribution <- function(sig,
                                bins = 40L,
                                observed_color = "#D55E00") {
    if (!inherits(sig, "decay_significance")) {
        stop("`sig` must be a `decay_significance` object from ",
            "test_decay_significance().")
    }
    .validate_positive_integer(bins, "bins")

    null_df <- data.frame(F_stat = sig$null_F)
    observed <- sig$anchor_test$F_stat

    subtitle <- sprintf(
        "Anchor: %s | n_random = %d successful fits | empirical p = %s",
        format(sig$anchor_position),
        length(sig$null_F),
        format.pval(sig$empirical_p, digits = 3)
    )

    ggplot2::ggplot(null_df, ggplot2::aes(x = F_stat)) +
        ggplot2::geom_histogram(bins = bins, fill = "grey70",
            color = "grey30") +
        ggplot2::geom_vline(xintercept = observed,
            color = observed_color, linewidth = 1.1) +
        ggplot2::labs(
            title    = "Null distribution of interaction F-statistics",
            subtitle = subtitle,
            x        = "Kenward-Roger F for distance x group",
            y        = "Count"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold")
        )
}

#' Plot coverage profile for a single sample
#'
#' Draws one sample's per-bin coverage against distance from the anchor,
#' with a dashed dark grey least-squares line as a visual guide. This is a
#' per-sample companion to [plot_coverage_decay()], useful for checking
#' that the decay in an individual sample looks like what the group
#' trend suggests.
#'
#' @details
#'
#' The fitted line is an ordinary least-squares fit of the plotted
#' response on distance from the anchor, using all of the sample's bins.
#' When `window` is supplied it only restricts the visible x-axis (via
#' `coord_cartesian`); the fit still uses bins outside the window, as in
#' [plot_coverage_decay()]. With `log_scale = TRUE` the response is the
#' `log_norm_cov` column, the same quantity the formal tests model, so
#' the line is the log-linear decay fit. With `log_scale = FALSE` the
#' response is `norm_cov` and the line is a linear fit on the raw scale,
#' which is not the model used for inference.
#'
#' The coverage line is colored by the sample's group, using the same
#' palette positions as the other plot functions, so a group has the same
#' color in every plot. The dashed fit line is always dark grey
#' (`"#4a4a4a"`) so it stays visible against any group color.
#'
#' If `combined` already contains a `distance_from_anchor` column, it is
#' used as-is and `anchor_position` is not evaluated. Otherwise
#' `anchor_position` is required and distances are computed via
#' [add_distance_from_anchor()].
#'
#' @param combined A data.frame or tibble from [load_all_samples()].
#'   Must contain columns `sample`, `group`, `bin_mid`, and `norm_cov`,
#'   plus `log_norm_cov` when `log_scale = TRUE`.
#' @param sample Character scalar. Name of the sample to plot. Must be
#'   present in `combined$sample` and belong to a single group.
#' @param anchor_position Non-negative numeric scalar. Reference
#'   coordinate for distance calculation. No default. Not needed when
#'   `combined` already contains `distance_from_anchor`.
#' @param window Optional positive numeric scalar. If supplied,
#'   restricts the x-axis to `[0, window]`. Default `NULL` shows the
#'   whole genome.
#' @param log_scale Logical. If `TRUE` (default), plot `log_norm_cov`; if
#'   `FALSE`, plot `norm_cov`.
#' @param palette Character scalar. `"okabe_ito"` (default) uses the
#'   Okabe-Ito colorblind-safe palette. A Brewer palette name (e.g.
#'   `"Set1"`) uses [ggplot2::scale_color_brewer()] instead. Any other
#'   string falls back to Okabe-Ito with a warning.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_coverage_decay()], [fit_decay_slopes()].
#'
#' @examples
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
#' plot_sample(combined, "A1", anchor_position = 1000, window = 1000)
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_smooth coord_cartesian labs theme_minimal theme element_text
#' @export
plot_sample <- function(combined,
                        sample,
                        anchor_position,
                        window = NULL,
                        log_scale = TRUE,
                        palette = "okabe_ito") {
    if (!is.data.frame(combined)) {
        stop("`combined` must be a data.frame or tibble.")
    }
    if (!is.logical(log_scale) || length(log_scale) != 1L ||
        is.na(log_scale)) {
        stop("`log_scale` must be TRUE or FALSE.")
    }
    required <- c("sample", "group", "bin_mid", "norm_cov",
        if (log_scale) "log_norm_cov")
    missing_cols <- setdiff(required, names(combined))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "`combined` is missing required column(s): %s",
            paste(missing_cols, collapse = ", ")
        ))
    }
    if (!is.character(sample) || length(sample) != 1L || is.na(sample)) {
        stop("`sample` must be a single character string.")
    }
    if (!sample %in% combined$sample) {
        stop(sprintf("`sample` '%s' not found in `combined`.", sample))
    }
    if (!is.null(window) &&
        (!is.numeric(window) || length(window) != 1L || window <= 0)) {
        stop("`window` must be a single positive number.")
    }
    if (!"distance_from_anchor" %in% names(combined)) {
        if (missing(anchor_position) || is.null(anchor_position)) {
            stop("Provide `anchor_position` or a `distance_from_anchor` ",
                "column in `combined`.")
        }
        combined <- add_distance_from_anchor(combined, anchor_position)
    }

    ## Factor over the full table so the group keeps its palette position
    ## from the other plots even though only one sample is drawn.
    combined$group <- factor(combined$group)
    sdat <- combined[combined$sample == sample, ]
    grp <- unique(as.character(sdat$group))
    if (length(grp) > 1L) {
        stop(sprintf(
            "`sample` '%s' belongs to multiple groups (%s); ",
            sample, paste(grp, collapse = ", ")
        ), "each sample must belong to exactly one group.")
    }

    y_col <- if (log_scale) "log_norm_cov" else "norm_cov"
    y_lab <- if (log_scale) "Log normalized coverage" else "Normalized coverage"
    subtitle <- if (log_scale) {
        "log-normalized coverage, log-linear fit shown"
    } else {
        "normalized coverage, linear fit shown"
    }

    p <- ggplot2::ggplot(
        sdat,
        ggplot2::aes(x = distance_from_anchor, y = .data[[y_col]],
            color = group)
    ) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_smooth(
            method = "lm", formula = y ~ x, linetype = "dashed", se = FALSE,
            color = "#4a4a4a"
        ) +
        .palette_scale(palette, nlevels(combined$group)) +
        ggplot2::labs(
            title    = paste0("Coverage profile: ", sample, " (", grp, ")"),
            subtitle = subtitle,
            x        = "Distance from anchor (bp)",
            y        = y_lab
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            plot.title      = ggplot2::element_text(face = "bold"),
            legend.position = "none"
        )
    if (!is.null(window)) {
        p <- p + ggplot2::coord_cartesian(xlim = c(0, window))
    }
    p
}
