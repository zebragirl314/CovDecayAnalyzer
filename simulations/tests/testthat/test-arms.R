## Contract checks for the A2-3 arm variants.
test_that("permutation Dunn equals omnibus permutation p at G = 2", {
    set.seed(11)
    g <- factor(rep(1:2, each = 4))
    x <- rnorm(8)
    r <- kw_perm_p(x, g)
    expect_equal(r$method, "exact_enum")
    expect_equal(unname(r$dunn_perm), r$p)
})

test_that("permutation Dunn runs by enumeration and by random draws", {
    set.seed(12)
    g3 <- factor(rep(1:3, each = 3))
    r3 <- kw_perm_p(rnorm(9), g3)
    expect_equal(r3$method, "exact_enum")
    expect_length(r3$dunn_perm, 3L)
    expect_true(all(is.finite(r3$dunn_perm)))
    g5 <- factor(rep(1:3, each = 5))
    r5 <- kw_perm_p(rnorm(15), g5, seed = 3)
    expect_equal(r5$method, "random_perm")
    expect_equal(kw_perm_p(rnorm(15), g5, seed = 3)$method, "random_perm")
})

test_that("pooled two-stage equals Welch arm at G = 3 and the F test at G = 2", {
    sc3 <- make_scenario(G = 3L, n = 4L, Delta_std = 1)
    d3 <- generate_replicate(sc3, 1L)
    expect_equal(analyze_two_stage_pooled(d3)$p_value,
        analyze_two_stage(d3)$p_value)
    sc2 <- make_scenario(G = 2L, n = 3L, Delta_std = 1)
    d2 <- generate_replicate(sc2, 1L)
    sl <- as.data.frame(CovDecayAnalyzer::fit_decay_slopes(d2))
    f <- anova(lm(slope ~ factor(group), data = sl))[["Pr(>F)"]][1L]
    expect_equal(analyze_two_stage_pooled(d2)$p_value, f)
})
