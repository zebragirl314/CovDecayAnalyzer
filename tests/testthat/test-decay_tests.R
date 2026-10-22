## Tests for decay_tests.R.
## These fit real lme4 models and can take a couple of seconds each. They
## are gated on lme4/lmerTest/emmeans being available so that the test suite
## remains runnable on stripped-down installations.

skip_if_not_installed("lme4")
skip_if_not_installed("lmerTest")
skip_if_not_installed("emmeans")


.make_fake_combined_lmm <- function(n_per_group = 5L, n_bins = 60L, seed = 1L,
                                    slope_g1 = -1e-5, slope_g2 = -3e-5,
                                    sample_re_sd = 5e-6,
                                    noise_sd = 0.12) {
    set.seed(seed)
    bin_mid <- seq(500, by = 1000, length.out = n_bins)
    rows <- list()
    idx <- 1L
    for (grp in c("g1", "g2")) {
        base_slope <- if (grp == "g1") slope_g1 else slope_g2
        for (i in seq_len(n_per_group)) {
            samp_slope <- base_slope + rnorm(1, sd = sample_re_sd)
            samp_intercept <- rnorm(1, sd = 0.1)
            rows[[idx]] <- tibble::tibble(
                sample = paste0(grp, "_", i),
                group  = grp,
                bin_mid = bin_mid,
                log_norm_cov = samp_intercept + samp_slope * bin_mid +
                    rnorm(n_bins, sd = noise_sd)
            )
            idx <- idx + 1L
        }
    }
    do.call(rbind, rows)
}


test_that("test_decay_slopes returns a decay_test object with expected fields", {
    fake <- .make_fake_combined_lmm()
    res <- test_decay_slopes(fake, anchor_position = 0)

    expect_s3_class(res, "decay_test")
    expect_true(all(c("model", "anova_table", "F_stat", "df1", "df2",
        "p_value", "singular", "conv_code", "group_slopes", "pairwise",
        "n_samples", "n_bins")
        %in% names(res)))
    expect_true(is.finite(res$F_stat))
    expect_true(is.finite(res$p_value))
    expect_gte(res$p_value, 0)
    expect_lte(res$p_value, 1)
})


test_that("test_decay_slopes detects a real group difference", {
    fake <- .make_fake_combined_lmm(
        n_per_group = 5L, n_bins = 60L,
        slope_g1 = -1e-5, slope_g2 = -5e-5,
        sample_re_sd = 3e-6, noise_sd = 0.10
    )
    res <- test_decay_slopes(fake, anchor_position = 0)
    expect_lt(res$p_value, 0.05)
})


test_that("test_decay_slopes does not flag equivalent groups", {
    fake <- .make_fake_combined_lmm(
        n_per_group = 5L, n_bins = 60L,
        slope_g1 = -1e-5, slope_g2 = -1e-5,
        sample_re_sd = 3e-6, noise_sd = 0.10
    )
    res <- test_decay_slopes(fake, anchor_position = 0)
    ## Should not consistently reject the null. Loosely: p > 0.01.
    expect_gt(res$p_value, 0.01)
})


test_that("test_decay_slopes requires at least two groups", {
    fake <- .make_fake_combined_lmm()
    fake$group <- "g1"
    expect_error(
        test_decay_slopes(fake, anchor_position = 0),
        "at least two groups"
    )
})


test_that("test_decay_slopes errors when a group has fewer than two samples", {
    fake <- .make_fake_combined_lmm(n_per_group = 3L)
    ## Keep only one g2 sample.
    one_g2 <- fake[fake$group == "g1" | fake$sample == "g2_1", ]
    expect_error(
        test_decay_slopes(one_g2, anchor_position = 0),
        "Group 'g2' has only 1 sample; at least 2 samples per group are required."
    )

    ## Several offending groups are all named, with their counts.
    fake3 <- .make_fake_combined_lmm(n_per_group = 3L)
    extra <- fake3[fake3$sample == "g1_1", ]
    extra$sample <- "g3_1"
    extra$group <- "g3"
    thin <- rbind(fake3[fake3$group == "g1" | fake3$sample == "g2_1", ],
        extra)
    expect_error(
        test_decay_slopes(thin, anchor_position = 0),
        "Group 'g2' has only 1 sample; Group 'g3' has only 1 sample"
    )
})


test_that("test_decay_slopes counts samples, not bins, per group", {
    ## Two samples per group with many bins each must be accepted.
    fake <- .make_fake_combined_lmm(n_per_group = 2L)
    res <- test_decay_slopes(fake, anchor_position = 0)
    expect_s3_class(res, "decay_test")
})


test_that("test_decay_slopes reports singularity and convergence diagnostics", {
    fake <- .make_fake_combined_lmm()
    res <- test_decay_slopes(fake, anchor_position = 0)
    expect_true(is.logical(res$singular) && length(res$singular) == 1L)
    expect_false(is.na(res$singular))
    expect_true(is.null(res$conv_code) || is.integer(res$conv_code))
})


test_that("print.decay_test shows a diagnostic only for singular or non-converged fits", {
    fake <- .make_fake_combined_lmm()
    res <- test_decay_slopes(fake, anchor_position = 0)

    clean <- res
    clean$singular <- FALSE
    clean$conv_code <- 0L
    out_clean <- capture.output(print(clean))
    expect_false(any(grepl("Warning:", out_clean)))

    sing <- res
    sing$singular <- TRUE
    sing$conv_code <- 0L
    out_sing <- capture.output(print(sing))
    expect_true(any(grepl(
        "Warning: model fit was singular; results should be interpreted with caution.",
        out_sing, fixed = TRUE
    )))

    nonconv <- res
    nonconv$singular <- FALSE
    nonconv$conv_code <- 1L
    out_nc <- capture.output(print(nonconv))
    expect_true(any(grepl("convergence code 1", out_nc)))
})
