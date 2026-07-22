## Tests for null_distribution.R.

test_that("empirical_pvalue is bounded and monotonic", {
    null_dist <- rf(500, df1 = 3, df2 = 40)

    p_low  <- empirical_pvalue(observed = 0.1, null_dist = null_dist)
    p_mid  <- empirical_pvalue(observed = 2,   null_dist = null_dist)
    p_high <- empirical_pvalue(observed = 100, null_dist = null_dist)

    expect_gte(p_low, 1 / (length(null_dist) + 1))
    expect_lte(p_low, 1)
    expect_lte(p_high, p_mid)
    expect_lte(p_mid, p_low)
    ## Rare observation gets the minimum p.
    expect_equal(p_high, 1 / (length(null_dist) + 1))
})


test_that("empirical_pvalue handles NAs in null", {
    null_dist <- c(1, 2, 3, NA, NA_real_)
    p <- empirical_pvalue(observed = 2.5, null_dist = null_dist)
    expect_equal(p, (1 + 1) / (3 + 1))
})


test_that("empirical_pvalue rejects empty finite null", {
    expect_error(
        empirical_pvalue(observed = 1, null_dist = c(NA_real_, NA_real_)),
        "no finite values"
    )
})


## Cheap smoke test that null_distribution_F wires together correctly.
## We use a very small n_random with skip_on_cran-worthy timing.
skip_if_not_installed("lme4")
skip_if_not_installed("lmerTest")
skip_if_not_installed("emmeans")


test_that("null_distribution_F returns a numeric vector under simulated data", {
    set.seed(1)
    fake <- do.call(rbind, lapply(seq_len(8), function(i) {
        grp <- if (i <= 4) "g1" else "g2"
        tibble::tibble(
            sample = paste0("s", i), group = grp,
            bin_mid = seq(500, by = 1000, length.out = 40),
            log_norm_cov = rnorm(40, sd = 0.2)
        )
    }))

    null_F <- null_distribution_F(
        fake, genome_length = 40000L, n_random = 5L, seed = 42
    )
    expect_type(null_F, "double")
    expect_lte(length(null_F), 5L)
    expect_true(all(is.finite(null_F)))
})
