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


test_that("null_distribution_F validates its arguments", {
    fake <- .make_fake_null_combined()
    expect_error(
        null_distribution_F("not a df", genome_length = 2000L),
        "data.frame or tibble"
    )
    expect_error(
        null_distribution_F(fake, genome_length = 0),
        "genome_length"
    )
    expect_error(
        null_distribution_F(fake, genome_length = 2000.5),
        "genome_length"
    )
    expect_error(
        null_distribution_F(fake, genome_length = 2000L, n_random = 0),
        "n_random"
    )
    expect_error(
        null_distribution_F(fake, genome_length = 2000L, n_random = 2.5),
        "n_random"
    )
    bad_exclude <- list("a", 1, c(1, 2, 3), c(5, 5), c(10, 5), c(1, NA),
        c(1, Inf))
    for (ex in bad_exclude) {
        expect_error(
            null_distribution_F(fake, genome_length = 2000L, n_random = 2L,
                exclude = ex),
            "exclude"
        )
    }
    expect_error(
        null_distribution_F(fake, genome_length = 2000L, n_random = 2L,
            seed = "a"),
        "seed"
    )
    expect_error(
        null_distribution_F(fake, genome_length = 2000L, n_random = 2L,
            seed = c(1, 2)),
        "seed"
    )
})


test_that("null_distribution_F reports failed fits in the n_failed attribute", {
    fake <- .make_fake_null_combined()

    ## Normal data: no failures, attribute present and zero.
    ok <- null_distribution_F(fake, genome_length = 2000L, n_random = 5L,
        seed = 1)
    expect_identical(attr(ok, "n_failed"), 0L)

    ## A single group makes every refit error, so every anchor fails.
    one_group <- fake
    one_group$group <- "g1"
    res <- null_distribution_F(one_group, genome_length = 2000L,
        n_random = 6L, seed = 1)
    expect_type(res, "double")
    expect_length(res, 0L)
    expect_identical(attr(res, "n_failed"), 6L)
})


test_that("null_distribution_F counts kept plus failed fits as the anchors drawn", {
    fake <- .make_fake_null_combined()
    res <- null_distribution_F(fake, genome_length = 2000L, n_random = 7L,
        seed = 3)
    expect_equal(length(res) + attr(res, "n_failed"), 7L)
})


test_that("test_decay_significance runs end to end and records failures", {
    fake <- .make_fake_null_combined()
    sig <- test_decay_significance(
        fake, anchor_position = 1000, genome_length = 2000L,
        n_random = 5L, seed = 11
    )
    expect_s3_class(sig, "decay_significance")
    expect_s3_class(sig$anchor_test, "decay_test")
    expect_true(is.numeric(sig$null_F))
    expect_null(attr(sig$null_F, "n_failed"))
    expect_identical(sig$n_null_failed, 0L)
    expect_gte(sig$empirical_p, 1 / (length(sig$null_F) + 1))
    expect_lte(sig$empirical_p, 1)
    expect_equal(sig$n_random, 5L)
    expect_output(print(sig), "5 of 5 random anchors succeeded")
})


test_that("test_decay_significance validates its arguments", {
    fake <- .make_fake_null_combined()
    expect_error(
        test_decay_significance(fake, anchor_position = -1,
            genome_length = 2000L, n_random = 2L),
        "anchor_position"
    )
    expect_error(
        test_decay_significance(fake, anchor_position = 3000,
            genome_length = 2000L, n_random = 2L),
        "exceeds"
    )
    expect_error(
        test_decay_significance(fake, anchor_position = 1000,
            genome_length = 2000L, n_random = 2L, exclude_radius = -5),
        "exclude_radius"
    )
    expect_error(
        test_decay_significance(fake, anchor_position = 1000,
            genome_length = 2000L, n_random = 2L, exclude_radius = c(1, 2)),
        "exclude_radius"
    )
})


test_that("test_decay_significance defaults exclude_radius to window", {
    fake <- .make_fake_null_combined()
    ## window = 990 around anchor 1000 excludes [10, 1990], leaving the 19
    ## positions 1:9 and 1991:2000. Asking for 25 must truncate to 19.
    expect_warning(
        sig <- test_decay_significance(
            fake, anchor_position = 1000, genome_length = 2000L,
            n_random = 25L, window = 990, seed = 5
        ),
        "25.*19"
    )
    expect_equal(length(sig$null_F) + sig$n_null_failed, 19L)

    ## Passing exclude_radius = 0 recovers the old no-exclusion behavior,
    ## so all 25 anchors are available and no truncation warning occurs.
    expect_no_warning(
        sig0 <- test_decay_significance(
            fake, anchor_position = 1000, genome_length = 2000L,
            n_random = 25L, window = 990, exclude_radius = 0, seed = 5
        )
    )
    expect_equal(length(sig0$null_F) + sig0$n_null_failed, 25L)
})


.make_printable_sig <- function(window = NULL, n_failed = 2L) {
    structure(
        list(
            anchor_test     = list(F_stat = 8.5, p_value = 0.002),
            null_F          = rf(18, df1 = 3, df2 = 40),
            n_null_failed   = n_failed,
            empirical_p     = 0.0526,
            anchor_position = 93000,
            genome_length   = 238080,
            n_random        = 20L,
            window          = window
        ),
        class = "decay_significance"
    )
}


test_that("print.decay_significance handles window NULL and window set", {
    set.seed(1)
    out_null <- capture.output(res <- print(.make_printable_sig()))
    expect_s3_class(res, "decay_significance")
    expect_true(any(grepl("full genome", out_null)))
    expect_true(any(grepl(
        "18 of 20 random anchors succeeded \\(2 dropped due to fit failure\\)",
        out_null
    )))

    out_win <- capture.output(print(.make_printable_sig(window = 50000)))
    expect_true(any(grepl("\\+/- 50000 bp", out_win)))
    expect_false(any(grepl("full genome", out_win)))
})


test_that("print.decay_significance falls back when n_null_failed is absent", {
    set.seed(1)
    sig <- .make_printable_sig()
    sig$n_null_failed <- NULL
    out <- capture.output(print(sig))
    expect_true(any(grepl("18 of 20 requested", out)))
})
