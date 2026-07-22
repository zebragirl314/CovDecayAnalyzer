## Tests for decay_slopes.R.

.make_fake_combined <- function(n_samples = 8L, n_bins = 50L, seed = 1L,
                                slope_g1 = -1e-5, slope_g2 = -3e-5,
                                noise_sd = 0.15) {
    set.seed(seed)
    bin_mid <- seq(500, by = 1000, length.out = n_bins)
    do.call(rbind, lapply(seq_len(n_samples), function(i) {
        grp <- if (i <= n_samples / 2) "g1" else "g2"
        rate <- if (grp == "g1") slope_g1 else slope_g2
        tibble::tibble(
            sample = paste0("s", i),
            group  = grp,
            bin_mid = bin_mid,
            log_norm_cov = rate * bin_mid + rnorm(n_bins, sd = noise_sd)
        )
    }))
}


test_that("add_distance_from_anchor computes absolute distance", {
    fake <- .make_fake_combined(n_samples = 2L, n_bins = 5L)
    out <- add_distance_from_anchor(fake, anchor_position = 2500)

    expect_true("distance_from_anchor" %in% names(out))
    expect_equal(out$distance_from_anchor, abs(out$bin_mid - 2500))
})


test_that("add_distance_from_anchor rejects negative anchor", {
    fake <- .make_fake_combined(n_samples = 2L, n_bins = 5L)
    expect_error(
        add_distance_from_anchor(fake, anchor_position = -100),
        "non-negative"
    )
})


test_that("fit_decay_slopes returns one row per sample", {
    fake <- .make_fake_combined()
    slopes <- fit_decay_slopes(fake, anchor_position = 0)

    expect_s3_class(slopes, "tbl_df")
    expect_equal(nrow(slopes), length(unique(fake$sample)))
    expect_true(all(c("sample", "group", "slope", "se_slope",
        "df", "n_bins") %in% names(slopes)))
})


test_that("fit_decay_slopes recovers the simulated slopes", {
    fake <- .make_fake_combined(n_samples = 20L, n_bins = 100L,
        slope_g1 = -1e-5, slope_g2 = -3e-5, noise_sd = 0.05)
    slopes <- fit_decay_slopes(fake, anchor_position = 0)

    g1_mean <- mean(slopes$slope[slopes$group == "g1"])
    g2_mean <- mean(slopes$slope[slopes$group == "g2"])

    expect_lt(abs(g1_mean - (-1e-5)), 5e-6)
    expect_lt(abs(g2_mean - (-3e-5)), 5e-6)
    expect_lt(g2_mean, g1_mean)
})


test_that("fit_decay_slopes windows correctly", {
    fake <- .make_fake_combined(n_samples = 4L, n_bins = 100L)
    ## use anchor 0 so distance == bin_mid; then window = 20000 keeps
    ## exactly the first 20 bins per sample.
    slopes <- fit_decay_slopes(fake, anchor_position = 0, window = 20000)
    expect_true(all(slopes$n_bins == 20L))
})
