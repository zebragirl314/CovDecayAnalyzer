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
