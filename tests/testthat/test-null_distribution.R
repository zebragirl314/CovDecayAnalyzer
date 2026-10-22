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


## Anchor-position sampling. null_distribution_F returns only F values, so
## these tests check the number of fits retained and that warnings fire (or
## do not) as appropriate. Small data keeps each refit cheap.
.make_fake_null_combined <- function(n_samples = 8L, n_bins = 20L,
                                     bin_width = 100) {
    set.seed(1)
    do.call(rbind, lapply(seq_len(n_samples), function(i) {
        tibble::tibble(
            sample = paste0("s", i),
            group = if (i <= n_samples / 2) "g1" else "g2",
            bin_mid = seq(bin_width / 2, by = bin_width,
                length.out = n_bins),
            log_norm_cov = rnorm(n_bins, sd = 0.2)
        )
    }))
}


test_that("null_distribution_F truncates with a warning when exclude leaves too few positions", {
    fake <- .make_fake_null_combined()
    ## exclude = c(1, 1900) leaves positions 1901..2000, i.e. 100 of them.
    expect_warning(
        res <- null_distribution_F(
            fake, genome_length = 2000L, n_random = 200L,
            exclude = c(1, 1900), seed = 1
        ),
        "200.*100"
    )
    expect_type(res, "double")
    expect_lte(length(res), 100L)
})


test_that("null_distribution_F truncates with a warning when n_random exceeds genome_length", {
    fake <- .make_fake_null_combined()
    expect_warning(
        res <- null_distribution_F(
            fake, genome_length = 100L, n_random = 500L, seed = 1
        ),
        "500.*100"
    )
    expect_type(res, "double")
    expect_lte(length(res), 100L)
    ## Every position is used, so nearly all fits should be retained.
    expect_gt(length(res), 0L)
})


test_that("null_distribution_F returns n_random fits with no warning in the normal case", {
    fake <- .make_fake_null_combined()
    expect_no_warning(
        res <- null_distribution_F(
            fake, genome_length = 2000L, n_random = 10L, seed = 1
        )
    )
    expect_length(res, 10L)

    expect_no_warning(
        res_ex <- null_distribution_F(
            fake, genome_length = 2000L, n_random = 10L,
            exclude = c(500, 1500), seed = 1
        )
    )
    expect_length(res_ex, 10L)
})


test_that("null_distribution_F is reproducible under a fixed seed", {
    fake <- .make_fake_null_combined()
    run <- function(seed, exclude = NULL) {
        null_distribution_F(
            fake, genome_length = 2000L, n_random = 8L,
            exclude = exclude, seed = seed
        )
    }
    expect_equal(run(7), run(7))
    expect_equal(run(7, c(500, 1500)), run(7, c(500, 1500)))
    ## Different seeds should sample different anchors, hence different F.
    expect_false(isTRUE(all.equal(run(7), run(8))))
})
