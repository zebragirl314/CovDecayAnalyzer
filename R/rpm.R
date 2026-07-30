## =========================================================================
## Viral reads-per-million (RPM) computation and visualization.
##
## The RPM analysis complements the coverage-decay analysis by
## summarizing overall viral read yield per sample as a scale-normalized
## quantity that can be compared across sequencing runs of different
## depths.
## =========================================================================


#' Compute viral reads-per-million (RPM)
#'
#' Vectorized calculation of viral RPM given per-sample viral read
#' counts and per-sample human read counts. RPM is defined as
#' `(viral_reads / human_reads) * 1e6`.
#'
#' @details
#'
#' RPM here normalizes viral read yield against sample-specific
#' sequencing depth. It answers "for every million reads that mapped
#' to the human host, how many mapped to the virus?" This is more
#' informative than raw viral read counts when comparing samples with
#' different total sequencing depths.
#'
#' In the Bosco lab pipeline, `viral_reads` typically comes from
#' breseq's `output/summary.json` field
#' `.reads.total_reads`, and `human_reads` typically comes from
#' the concordant read count in bowtie2's alignment log
#' (`Illumina_GRCh38.txt`). Both quantities are extracted per sample
#' by the accompanying shell scripts and combined into a TSV that
#' can be read into R with `read.delim()` or [utils::read.table()].
#'
#' If any element of `human_reads` is zero, the corresponding RPM is
#' returned as `NA_real_` and a warning is issued.
#'
#' @param viral_reads Numeric vector. Number of viral reads per sample.
#'   Must be non-negative and the same length as `human_reads`.
#' @param human_reads Numeric vector. Number of human reads per sample.
#'   Must be non-negative and the same length as `viral_reads`.
#'
#' @return Numeric vector of RPM values, one per input sample. Values
#'   are `NA_real_` for samples with zero human reads.
#'
#' @seealso [plot_rpm_by_group()] for the standard summary
#'   visualization.
#'
#' @examples
#' compute_viral_rpm(
#'     viral_reads = c(12345, 67890, 45000),
#'     human_reads = c(1e7, 1.2e7, 9e6)
#' )
#'
#' @export
compute_viral_rpm <- function(viral_reads, human_reads) {
    if (!is.numeric(viral_reads) || !is.numeric(human_reads)) {
        stop("`viral_reads` and `human_reads` must be numeric vectors.")
    }
    if (length(viral_reads) != length(human_reads)) {
        stop("`viral_reads` and `human_reads` must be the same length.")
    }
    if (any(viral_reads < 0, na.rm = TRUE) ||
        any(human_reads < 0, na.rm = TRUE)) {
        stop("`viral_reads` and `human_reads` must be non-negative.")
    }

    rpm <- (viral_reads / human_reads) * 1e6
    zero_human <- !is.na(human_reads) & human_reads == 0
    if (any(zero_human)) {
        warning(sprintf(
            "Zero human reads in %d sample(s); RPM set to NA.",
            sum(zero_human)
        ))
        rpm[zero_human] <- NA_real_
    }
    rpm
}


#' Barplot of viral RPM by treatment group
#'
#' Draws a barplot of mean viral RPM per treatment group with
#' error bars, optionally overlaying individual sample points.
#'
#' @details
#'
#' The default error bar is `mean +/- SD`, which is a plain summary of
#' within-group spread. Users who want a formal comparison should
#' first fit a suitable model (e.g. `lm(log(viral_rpm) ~ group)`) and
#' report those results in a separate table, since the barplot is a
#' visual summary, not a test.
#'
#' Individual sample points overlaid on the bars help readers judge
#' whether the group mean is dominated by a single high-RPM sample
#' (an important sanity check in small-n experiments). Set
#' `show_points = FALSE` for a cleaner presentation when n per group
#' is large.
#'
#' @param rpm_df A data.frame or tibble containing per-sample RPM
#'   values. Must contain the columns named by `group_col` and
#'   `rpm_col`. Typically produced by joining a sample metadata table
#'   with a vector of RPM values from [compute_viral_rpm()].
#' @param group_col Character scalar. Name of the column containing
#'   the treatment group label. Default `"group"`.
#' @param rpm_col Character scalar. Name of the column containing the
#'   viral RPM values. Default `"viral_rpm"`.
#' @param show_points Logical. If `TRUE` (default), overlay each
#'   sample as a jittered point on top of the group bar.
#' @param bar_fill Character scalar. Fill color for the bars. Default
#'   `"#00693E"` (Dartmouth green), matching the prior lab convention.
#' @param error_color Character scalar. Color of the error bar
#'   segments. Default `"orange"`, matching the prior lab convention.
#'
#' @return A `ggplot` object.
#'
#' @seealso [compute_viral_rpm()].
#'
#' @examples
#' set.seed(1)
#' rpm_df <- data.frame(
#'     sample     = paste0("s", 1:12),
#'     group      = rep(c("WT-GCV", "WT+GCV", "FV-GCV", "FV+GCV"),
#'                      each = 3),
#'     viral_rpm  = c(rnorm(3, 3000, 500), rnorm(3, 1500, 300),
#'                    rnorm(3, 2500, 400), rnorm(3, 1200, 250))
#' )
#' plot_rpm_by_group(rpm_df)
#'
#' @importFrom ggplot2 ggplot aes stat_summary geom_jitter mean_sdl labs theme_bw theme element_text
#' @export
plot_rpm_by_group <- function(rpm_df,
                            group_col = "group",
                            rpm_col   = "viral_rpm",
                            show_points = TRUE,
                            bar_fill    = "#00693E",
                            error_color = "orange") {
    if (!is.data.frame(rpm_df)) {
        stop("`rpm_df` must be a data.frame or tibble.")
    }
    .validate_string(group_col, "group_col")
    .validate_string(rpm_col, "rpm_col")
    missing_cols <- setdiff(c(group_col, rpm_col), names(rpm_df))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "`rpm_df` is missing required column(s): %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    rpm_df[[group_col]] <- factor(rpm_df[[group_col]])

    p <- ggplot2::ggplot(
        rpm_df,
        ggplot2::aes(x = .data[[group_col]], y = .data[[rpm_col]])
    ) +
        ggplot2::stat_summary(
            fun = mean, geom = "bar",
            fill = bar_fill, alpha = 0.8, width = 0.7
        ) +
        ggplot2::stat_summary(
            fun.data = ggplot2::mean_sdl,
            fun.args = list(mult = 1),
            geom = "errorbar", width = 0.2,
            colour = error_color, linewidth = 1
        ) +
        ggplot2::labs(
            y = "Viral RPM",
            x = NULL,
            title = "Viral reads per million (RPM) by treatment group"
        ) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.title  = ggplot2::element_text(face = "bold")
        )

    if (isTRUE(show_points)) {
        p <- p + ggplot2::geom_jitter(width = 0.15, size = 2,
            alpha = 0.75, colour = "black")
    }
    p
}
