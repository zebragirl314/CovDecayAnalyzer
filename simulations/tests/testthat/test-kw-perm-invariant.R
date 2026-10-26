## Invariant: the vectorized random-permutation path of kw_perm_p() and the
## preserved loop implementation (analyze_kw_perm_loop_internal) give the
## same null behaviour. G = 3, n = 5, Delta_std = 0, reference nuisance (the
## random-permutation path: 756,756 assignments exceed the enumeration cap),
## 500 datasets with distinct seeds.

.cp <- function(x, n) {
    ci <- stats::binom.test(x, n)$conf.int
    c(lower = ci[[1L]], upper = ci[[2L]])
}

## Overlap of two intervals as a fraction of the narrower one.
.overlap_frac <- function(a, b) {
    ov <- min(a[["upper"]], b[["upper"]]) - max(a[["lower"]], b[["lower"]])
    max(ov, 0) / min(a[["upper"]] - a[["lower"]], b[["upper"]] - b[["lower"]])
}

test_that("vectorized kw_perm matches the loop path on 500 null datasets", {
    sc <- make_scenario(G = 3L, n = 5L, Delta_std = 0)
    R <- 500L
    cores <- max(1L, min(10L, parallel::detectCores() - 2L))
    res <- parallel::mclapply(seq_len(R), function(r) {
        d <- generate_replicate(sc, r)
        new <- analyze_kw_perm(d)
        old <- analyze_kw_perm_loop_internal(d)
        list(new = new, old = old)
    }, mc.cores = cores)
    pick <- function(which, variant) {
        vapply(res, function(x) {
            r <- x[[which]]
            r$p_value[r$arm == variant]
        }, numeric(1L))
    }
    ## No fit failed and every fit used the random-permutation path.
    expect_true(all(vapply(res, function(x) all(x$new$status == "ok") &&
        all(x$old$status == "ok"), logical(1L))))
    expect_true(all(vapply(res, function(x) all(x$new$kw_perm_method ==
        "random_perm"), logical(1L))))

    for (v in c("kw_perm", "kw_perm_dunn_perm_ungated",
        "kw_perm_dunn_perm_gated")) {
        pn <- pick("new", v)
        po <- pick("old", v)
        xn <- sum(pn <= 0.05)
        xo <- sum(po <= 0.05)
        cn <- .cp(xn, R)
        co <- .cp(xo, R)
        ## The criterion from the task: CP intervals overlap by >= 50%.
        expect_gte(.overlap_frac(cn, co), 0.5)
        ## Stronger: the draws are the same stream, so the p-values are
        ## identical, not merely close.
        expect_equal(pn, po)
        message(sprintf("%s: new %d/%d (CP %.3f-%.3f), loop %d/%d (CP %.3f-%.3f), overlap %.2f, max |dp| %.2g",
            v, xn, R, cn[[1L]], cn[[2L]], xo, R, co[[1L]], co[[2L]],
            .overlap_frac(cn, co), max(abs(pn - po))))
    }
})

test_that("vectorized and loop paths agree on the statistics directly", {
    sc <- make_scenario(G = 4L, n = 4L, Delta_std = 1)
    d <- generate_replicate(sc, 11L)
    sl <- .sample_slopes(d)
    a <- kw_perm_p(sl$slope, sl$group, seed = 3L)
    b <- kw_perm_p_loop_internal(sl$slope, sl$group, seed = 3L)
    expect_identical(a, b)
    ## Exact-enumeration path is untouched.
    sl2 <- sl[1:6, ]
    sl2$group <- factor(rep(c("g1", "g2"), each = 3L))
    expect_equal(kw_perm_p(sl2$slope, sl2$group)$method, "exact_enum")
    expect_identical(kw_perm_p(sl2$slope, sl2$group),
        kw_perm_p_loop_internal(sl2$slope, sl2$group))
})
