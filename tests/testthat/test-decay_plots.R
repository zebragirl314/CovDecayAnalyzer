## Tests for decay_plots.R.
## We check that plots return ggplot objects and validate their inputs
## correctly. We do not render or snapshot plots because ggplot rendering
## is not reproducible across ggplot2 versions in a way that survives
## R CMD check.

.make_fake_combined_plots <- function() {
    set.seed(1)
    do.call(rbind, lapply(paste0("s", 1:6), function(sname) {
        grp <- if (as.integer(sub("s", "", sname)) <= 3) "g1" else "g2"
        rate <- if (grp == "g1") -1e-5 else -3e-5
        tibble::tibble(
            sample  = sname, group = grp,
            bin_mid = seq(500, 50000, by = 1000),
            norm_cov = exp(rate * seq(500, 50000, by = 1000)) +
                rnorm(50, sd = 0.05),
            log_norm_cov = rate * seq(500, 50000, by = 1000) +
                rnorm(50, sd = 0.05)
        )
    }))
}


test_that("plot_coverage_decay returns a ggplot object", {
    fake <- .make_fake_combined_plots()
    p <- plot_coverage_decay(fake, anchor_position = 0)
    expect_s3_class(p, "ggplot")
})


test_that("plot_coverage_decay errors on missing columns", {
    bad <- tibble::tibble(sample = "s1", group = "g1",
        bin_mid = 500, log_norm_cov = 0)  # no norm_cov
    expect_error(
        plot_coverage_decay(bad, anchor_position = 0),
        "norm_cov"
    )
})


test_that("plot_coverage_decay applies window via coord_cartesian", {
    fake <- .make_fake_combined_plots()
    p <- plot_coverage_decay(fake, anchor_position = 0, window = 20000)
    ## The layer should apply xlim = c(0, 20000) via coord_cartesian.
    expect_equal(p$coordinates$limits$x, c(0, 20000))
})


test_that("plot_decay_slopes returns a ggplot object", {
    fake <- .make_fake_combined_plots()
    slopes <- fit_decay_slopes(fake, anchor_position = 0)
    p <- plot_decay_slopes(slopes)
    expect_s3_class(p, "ggplot")
})


test_that("plot_decay_slopes errors when input lacks required columns", {
    expect_error(
        plot_decay_slopes(data.frame(other = 1:3)),
        "group"
    )
})


test_that("plot_null_distribution returns a ggplot object", {
    set.seed(1)
    fake_sig <- structure(
        list(
            anchor_test     = list(F_stat = 8.5, p_value = 0.002),
            null_F          = rf(200, df1 = 3, df2 = 40),
            empirical_p     = 0.006,
            anchor_position = 93000,
            genome_length   = 238080,
            n_random        = 200,
            window          = 50000
        ),
        class = "decay_significance"
    )
    p <- plot_null_distribution(fake_sig)
    expect_s3_class(p, "ggplot")
})


test_that("plot_null_distribution rejects non-decay_significance input", {
    expect_error(
        plot_null_distribution(list(F_stat = 1)),
        "decay_significance"
    )
})


.fake_sig_plots <- function() {
    set.seed(1)
    structure(
        list(
            anchor_test     = list(F_stat = 8.5, p_value = 0.002),
            null_F          = rf(200, df1 = 3, df2 = 40),
            empirical_p     = 0.006,
            anchor_position = 93000,
            genome_length   = 238080,
            n_random        = 200,
            window          = 50000
        ),
        class = "decay_significance"
    )
}


test_that("plot_coverage_decay defaults to the Okabe-Ito palette", {
    fake <- .make_fake_combined_plots()
    p <- plot_coverage_decay(fake, anchor_position = 0)
    ## Layer 2 is the loess smooth; one color per group, in level order.
    cols <- unique(ggplot2::layer_data(p, 2)$colour)
    expect_equal(sort(cols), sort(.okabe_ito[1:2]))
})


test_that("plot_coverage_decay still accepts Brewer palette names", {
    fake <- .make_fake_combined_plots()
    p <- plot_coverage_decay(fake, anchor_position = 0, palette = "Set1")
    cols <- unique(ggplot2::layer_data(p, 2)$colour)
    expect_equal(sort(cols), sort(c("#E41A1C", "#377EB8")))
})


test_that("an unknown palette name warns and falls back to Okabe-Ito", {
    fake <- .make_fake_combined_plots()
    expect_warning(
        p <- plot_coverage_decay(fake, anchor_position = 0,
            palette = "nonsense"),
        "Unknown palette"
    )
    cols <- unique(ggplot2::layer_data(p, 2)$colour)
    expect_equal(sort(cols), sort(.okabe_ito[1:2]))
})


test_that("the Okabe-Ito palette errors when there are too many groups", {
    fake <- .make_fake_combined_plots()
    fake$group <- paste0("g", rep(1:9, length.out = nrow(fake)))
    expect_error(
        plot_coverage_decay(fake, anchor_position = 0),
        "Okabe-Ito palette has 8 colors"
    )
})


test_that("plot_decay_slopes defaults to the Okabe-Ito palette", {
    fake <- .make_fake_combined_plots()
    slopes <- fit_decay_slopes(fake, anchor_position = 0)
    p <- plot_decay_slopes(slopes)
    expect_equal(sort(unique(ggplot2::layer_data(p, 1)$fill)),
        sort(.okabe_ito[1:2]))
})


test_that("plot_decay_slopes still accepts Brewer palette names", {
    fake <- .make_fake_combined_plots()
    slopes <- fit_decay_slopes(fake, anchor_position = 0)
    p <- plot_decay_slopes(slopes, palette = "Set1")
    expect_equal(sort(unique(ggplot2::layer_data(p, 1)$fill)),
        sort(c("#E41A1C", "#377EB8")))
})


test_that("plot_null_distribution observed-F line uses Okabe-Ito vermillion", {
    p <- plot_null_distribution(.fake_sig_plots())
    ## Layer 2 is the vertical line.
    expect_equal(ggplot2::layer_data(p, 2)$colour, "#D55E00")
})


test_that("plot_null_distribution observed_color can still be overridden", {
    p <- plot_null_distribution(.fake_sig_plots(), observed_color = "blue")
    expect_equal(ggplot2::layer_data(p, 2)$colour, "blue")
})


test_that("plot_sample returns a ggplot object", {
    fake <- .make_fake_combined_plots()
    p <- plot_sample(fake, "s1", anchor_position = 0)
    expect_s3_class(p, "ggplot")
})


test_that("plot_sample errors when a required column is missing", {
    fake <- .make_fake_combined_plots()
    expect_error(
        plot_sample(fake[, setdiff(names(fake), "norm_cov")], "s1",
            anchor_position = 0),
        "norm_cov"
    )
    ## log_norm_cov is only required when log_scale = TRUE.
    no_log <- fake[, setdiff(names(fake), "log_norm_cov")]
    expect_error(plot_sample(no_log, "s1", anchor_position = 0),
        "log_norm_cov")
    expect_s3_class(
        plot_sample(no_log, "s1", anchor_position = 0, log_scale = FALSE),
        "ggplot"
    )
})


test_that("plot_sample errors on an unknown sample", {
    fake <- .make_fake_combined_plots()
    expect_error(
        plot_sample(fake, "nope", anchor_position = 0),
        "not found"
    )
    expect_error(
        plot_sample(fake, c("s1", "s2"), anchor_position = 0),
        "single character string"
    )
})


test_that("plot_sample errors when a sample matches multiple groups", {
    fake <- .make_fake_combined_plots()
    fake$group[fake$sample == "s1"][1:3] <- "g2"
    expect_error(
        plot_sample(fake, "s1", anchor_position = 0),
        "multiple groups"
    )
})


test_that("plot_sample requires anchor_position unless distances are present", {
    fake <- .make_fake_combined_plots()
    expect_error(plot_sample(fake, "s1"), "anchor_position")
    with_dist <- add_distance_from_anchor(fake, 0)
    expect_s3_class(plot_sample(with_dist, "s1"), "ggplot")
})


test_that("plot_sample validates window and log_scale", {
    fake <- .make_fake_combined_plots()
    expect_error(
        plot_sample(fake, "s1", anchor_position = 0, window = -1),
        "window"
    )
    expect_error(
        plot_sample(fake, "s1", anchor_position = 0, log_scale = NA),
        "log_scale"
    )
})


test_that("plot_sample respects window and defaults to the whole genome", {
    fake <- .make_fake_combined_plots()
    p <- plot_sample(fake, "s1", anchor_position = 0, window = 20000)
    expect_equal(p$coordinates$limits$x, c(0, 20000))
    p_all <- plot_sample(fake, "s1", anchor_position = 0)
    expect_null(p_all$coordinates$limits$x)
})


test_that("plot_sample uses distance_from_anchor as the x aesthetic", {
    fake <- .make_fake_combined_plots()
    p <- plot_sample(fake, "s1", anchor_position = 10000)
    expect_equal(all.vars(p$mapping$x), "distance_from_anchor")
    ## x values are distances, not genomic coordinates.
    expect_equal(sort(unique(ggplot2::layer_data(p, 1)$x)),
        sort(unique(abs(fake$bin_mid[fake$sample == "s1"] - 10000))))
})


test_that("plot_sample plots log_norm_cov directly when log_scale is TRUE", {
    fake <- .make_fake_combined_plots()
    p <- plot_sample(fake, "s1", anchor_position = 0)
    expect_equal(ggplot2::layer_data(p, 1)$y,
        fake$log_norm_cov[fake$sample == "s1"][order(
            fake$bin_mid[fake$sample == "s1"])])
    p_raw <- plot_sample(fake, "s1", anchor_position = 0, log_scale = FALSE)
    expect_equal(ggplot2::layer_data(p_raw, 1)$y,
        fake$norm_cov[fake$sample == "s1"][order(
            fake$bin_mid[fake$sample == "s1"])])
})


test_that("plot_sample subtitle describes what is drawn", {
    fake <- .make_fake_combined_plots()
    expect_match(
        plot_sample(fake, "s1", anchor_position = 0)$labels$subtitle,
        "log-linear fit shown"
    )
    expect_match(
        plot_sample(fake, "s1", anchor_position = 0,
            log_scale = FALSE)$labels$subtitle,
        "linear fit shown"
    )
})


test_that("plot_sample colors by group using the same palette positions", {
    fake <- .make_fake_combined_plots()
    ## s1 is in g1 (first level), s6 in g2 (second level).
    p1 <- plot_sample(fake, "s1", anchor_position = 0)
    p6 <- plot_sample(fake, "s6", anchor_position = 0)
    expect_equal(unique(ggplot2::layer_data(p1, 1)$colour), .okabe_ito[1])
    expect_equal(unique(ggplot2::layer_data(p6, 1)$colour), .okabe_ito[2])
    p6_brewer <- plot_sample(fake, "s6", anchor_position = 0,
        palette = "Set1")
    expect_equal(unique(ggplot2::layer_data(p6_brewer, 1)$colour), "#377EB8")
})


test_that("plot_sample draws the fit line in dark grey", {
    fake <- .make_fake_combined_plots()
    p <- plot_sample(fake, "s1", anchor_position = 0)
    expect_equal(unique(ggplot2::layer_data(p, 2)$colour), "#4a4a4a")
})
