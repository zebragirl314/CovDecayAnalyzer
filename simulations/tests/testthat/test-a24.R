## A2-4 additions: Sobol design, Tier 2 grids, sensitivity hooks, equal-pair
## post-hoc p-values, and the between-stage check function.

test_that("Sobol points are dyadically stratified and evenly spread", {
    x <- rbind(0, sobol_points(511L, 7L))  # indices 0..511
    for (j in 1:7) {
        expect_equal(sort(floor(x[, j] * 512)), 0:511)
    }
    x4 <- sobol_points(400L, 7L)
    for (i in 1:6) for (j in (i + 1):7) {
        tb <- table(factor(floor(x4[, i] * 4), 0:3),
            factor(floor(x4[, j] * 4), 0:3))
        expect_lt(max(abs(tb - 25)), 8)
    }
    ## Known leading points of the Joe-Kuo sequence (dims 1 to 2).
    expect_equal(sobol_points(3L, 2L)[, 2L], c(0.5, 0.25, 0.75))
})

test_that("A2-4 grids have the specified sizes and unique scenarios", {
    g1 <- build_grid("a24typeI")
    expect_equal(c(nrow(g1), sum(g1$reps)), c(27L, 135000L))
    expect_equal(sum(g1$family == "typeI_partial"), 6L)
    g2 <- build_grid("a24power")
    expect_equal(sum(g2$family == "power_fig"), 63L)
    expect_equal(sum(g2$family == "power_surface"), 400L)
    expect_equal(sum(g2$reps[g2$family == "power_surface"]), 120000L)
    g3 <- build_grid("a24sens")
    expect_true(all(g3$reps == 1000L))
    expect_equal(sum(g3$family == "sens_ref"), 12L)
    for (g in list(g1, g2, g3)) {
        expect_false(anyDuplicated(g[c("scenario_id", "bin_size")]) > 0L)
    }
    ## Surface design stays inside the stated ranges.
    p <- a24_surface_points()
    expect_true(all(p$n %in% 2:12) && all(p$G %in% 2:4))
    expect_true(all(p$Delta_std > 0 & p$Delta_std < 2))
    expect_true(all(p$cv > 0.1 & p$cv < 0.5))
    expect_true(all(p$sigma_atom > 0.05 & p$sigma_atom < 0.5))
    expect_true(all(p$L >= 50000 & p$L <= 500000))
    expect_true(all(p$anchor_frac > 0.1 & p$anchor_frac < 0.9))
})

test_that("sensitivity cells pair with their reference through seed_id", {
    g3 <- build_grid("a24sens")
    ref <- g3[g3$family == "sens_ref", ]
    for (f in c("sens_S1a", "sens_S1b", "sens_S2", "sens_S3", "sens_S4")) {
        expect_true(all(g3$seed_id[g3$family == f] %in% ref$scenario_id))
    }
    ## Baseline scenarios keep their ID as seed.
    expect_equal(ref$seed_id, ref$scenario_id)
    sc <- scenario_from_row(g3[g3$family == "sens_S0", ][1L, ])
    expect_equal(sc$n, c(2L, 4L))
    sc2 <- scenario_from_row(g3[g3$family == "sens_S2", ][1L, ])
    expect_equal(sc2$sens$S2$frac, 0.005)
})

test_that("baseline replicates are unchanged by the A2-4 generator edits", {
    f <- file.path(.sim_root, "results", "explore", "tier1",
        "tier1-replicates.rds")
    skip_if_not(file.exists(f))
    t1 <- readRDS(f)
    sc <- make_scenario(G = 3L, n = 5L, Delta_std = 1)
    ref <- t1[t1$scenario_id == sc$scenario_id & t1$bin_size == 5000 &
        t1$arm == "lmm" & t1$replicate_id <= 2, ]
    ref <- ref[order(ref$replicate_id), ]
    new <- do.call(rbind, lapply(1:2, function(r) {
        analyze_lmm(generate_replicate(sc, r))
    }))
    expect_equal(new$p_value, ref$p_value)
})

test_that("sensitivity hooks do what they claim", {
    base <- make_scenario(G = 2L, n = 3L, Delta_std = 0.5)
    d0 <- generate_replicate(base, 7L)
    ## Neutral settings reproduce the baseline exactly.
    s1 <- make_scenario(G = 2L, n = 3L, Delta_std = 0.5,
        sens = list(S1a = list(phi = 0)))
    s3 <- make_scenario(G = 2L, n = 3L, Delta_std = 0.5,
        sens = list(S3 = list(sigma_ratio = 1, cv_ratio = 1)))
    expect_equal(generate_replicate(s1, 7L)$log_norm_cov, d0$log_norm_cov)
    expect_equal(generate_replicate(s3, 7L)$log_norm_cov, d0$log_norm_cov)
    ## AR(1): marginal SD preserved, lag-1 correlation about phi.
    e <- matrix(stats::rnorm(4000 * 5, sd = 0.15), 4000, 5)
    sc <- make_scenario(G = 2L, n = 2L, sens = list(S1a = list(phi = 0.6)))
    x <- sens_s1a_noise(e, sc)
    expect_equal(stats::sd(x), 0.15, tolerance = 0.03)
    expect_equal(mean(vapply(1:5, function(j) stats::cor(x[-1, j],
        x[-nrow(x), j]), 0)), 0.6, tolerance = 0.05)
    ## Dropouts: zeroed fraction near the target and exactly zero atoms.
    sc2 <- make_scenario(G = 2L, n = 5L,
        sens = list(S2 = list(frac = 0.05, run_kb = 5)))
    z <- vapply(1:30, function(i) {
        mean(with_seed(i, simulate_atoms(sc2))$cov == 0)
    }, 0)
    expect_equal(mean(z), 0.05, tolerance = 0.2)
    ## Floor: coverage far from the anchor stays near the floor.
    expect_equal(exp(sens_s1b_curve(1e7, 6e-6, make_scenario(G = 2L,
        sens = list(S1b = list(floor = 0.4))))), matrix(0.4), ignore_attr = TRUE)
    ## Heteroscedastic scales apply to the last group only.
    sc3 <- make_scenario(G = 3L, n = 2L,
        sens = list(S3 = list(sigma_ratio = 2, cv_ratio = 3)))
    tr <- attr(generate_replicate(sc3, 1L), "truth")
    expect_equal(unique(tr$sigma_i[tr$group == "g3"]), 0.3)
    expect_equal(unique(tr$cv_i[tr$group == "g3"]), 1)
    expect_equal(unique(tr$sigma_i[tr$group != "g3"]), 0.15)
    ## Bad parameters fail loudly.
    expect_error(sens_s1a_noise(e, make_scenario(G = 2L,
        sens = list(S1a = list(phi = 1.2)))))
})

test_that("unbalanced and end-anchor scenarios run through every arm", {
    d <- generate_replicate(make_scenario(G = 2L, n = c(2L, 4L),
        Delta_std = 1), 2L)
    expect_equal(length(unique(d$sample)), 6L)
    expect_equal(analyze_lmm(d)$status, "ok")
    expect_equal(analyze_two_stage_pooled(d)$status, "ok")
    expect_true(all(analyze_kw_perm(d)$status == "ok"))
    for (a in c(0, 1)) {
        d <- generate_replicate(make_scenario(G = 2L, n = 3L, anchor_frac = a), 2L)
        expect_equal(analyze_lmm(d)$status, "ok")
    }
})

test_that("equal-pair post-hoc p-values pick only truly equal pairs", {
    d <- generate_replicate(make_scenario(G = 3L, n = 4L, Delta_std = 1), 3L)
    kw <- analyze_kw_perm(d)
    dr <- kw[kw$arm == "kw_perm_dunn_ungated", ]
    expect_true(is.finite(dr$p_equal_min))
    ## Only g1-g2 is a truly equal pair when g3 carries the effect.
    sl <- .sample_slopes(d)
    pa <- .dunn_holm(sl$slope, sl$group)
    expect_equal(dr$p_equal_min, pa[["g1|g2"]], ignore_attr = TRUE)
    lm <- analyze_lmm(d)
    expect_true(lm$p_equal_min >= 0 && lm$p_equal_min <= 1)
    ## Global null: every pair is equal, so the minimum is over all pairs.
    d0 <- generate_replicate(make_scenario(G = 3L, n = 4L, Delta_std = 0), 3L)
    sl0 <- .sample_slopes(d0)
    expect_equal(.min_equal_p(d0, .dunn_holm(sl0$slope, sl0$group)),
        min(.dunn_holm(sl0$slope, sl0$group)))
})

.fake_stage <- function(n_rep = 400L, p_lmm = function(n) stats::runif(n),
                        delta = 0, pooled = NULL, fail = 0) {
    mk <- function(arm, p) data.frame(scenario_id = "sc1", bin_size = 5000L,
        replicate_id = seq_len(n_rep), arm = arm, test_variant = "omnibus",
        status = ifelse(seq_len(n_rep) <= fail, "boom", "ok"), p_value = p,
        singular = FALSE, family = "typeI_null", G = 2L, n = 3L,
        n_spec = "3", Delta_std = delta, sweep_factor = "reference",
        sweep_level = "reference", stringsAsFactors = FALSE)
    p <- p_lmm(n_rep)
    rbind(mk("lmm", p), mk("two_stage_pooled", if (is.null(pooled)) p else pooled))
}

test_that("check_stage_results flags what the prompt says to flag", {
    set.seed(11)
    ok <- check_stage_results("typeI", .fake_stage(2000L))
    expect_true(ok$continue)
    expect_equal(nrow(ok$hard_stop_flags), 0L)
    ## Type I far above 0.07: hard stop. At the relaxed bound it is still one.
    hi <- .fake_stage(2000L, p_lmm = function(n) stats::runif(n, 0, 0.5))
    expect_false(check_stage_results("typeI", hi)$continue)
    ## Type I far below 0.025: hard stop.
    lo <- .fake_stage(2000L, p_lmm = function(n) stats::runif(n, 0.1, 1))
    chk <- check_stage_results("typeI", lo)
    expect_true("lmm_typeI_ci_below_0.025" %in% chk$hard_stop_flags$flag)
    ## 0.075 passes the sensitivity bound (0.08) but not the main one.
    p75 <- .fake_stage(20000L, p_lmm = function(n) ifelse(seq_len(n) <= 1500L, 0.01, 0.5))
    expect_true("lmm_typeI_ci_above_upper" %in%
        check_stage_results("typeI", p75)$hard_stop_flags$flag)
    expect_false("lmm_typeI_ci_above_upper" %in%
        check_stage_results("sensitivities", p75)$hard_stop_flags$flag)
    ## Arm failures: 20% hard, 5% soft.
    expect_true("arm_failure_gt_10pct" %in%
        check_stage_results("typeI", .fake_stage(400L, fail = 80L))$hard_stop_flags$flag)
    sf <- check_stage_results("typeI", .fake_stage(400L, fail = 20L))
    expect_true("arm_failure_1_to_10pct" %in% sf$soft_stop_flags$flag)
    expect_true(sf$continue || nrow(sf$hard_stop_flags) > 0L)
    ## Pooled beating the LMM by more than 0.02 at a power cell: soft flag.
    pw <- .fake_stage(1000L, p_lmm = function(n) ifelse(seq_len(n) <= 500L, 0.01, 0.9),
        delta = 1, pooled = ifelse(seq_len(1000L) <= 600L, 0.01, 0.9))
    cp <- check_stage_results("power", pw)
    expect_true("pooled_two_stage_beats_lmm_gt_0.02" %in% cp$soft_stop_flags$flag)
    expect_true(cp$continue)
    ## Panic exits.
    rs <- data.frame(peak_worker_memory_mb = 4000, cache_size_gb = 12)
    pn <- check_stage_results("typeI", .fake_stage(400L), resource = rs,
        wall_hours = 11, budget_hours = 6.5, total_wall_hours = 25)
    expect_setequal(pn$panic_exit_flags$flag, c("worker_memory_gt_3GB",
        "cache_size_gt_10GB", "stage_wall_gt_150pct_budget", "total_wall_gt_cap"))
    expect_false(pn$continue)
})

test_that("run_batch can restrict arms (Sobol surface points are LMM-only)", {
    g <- build_grid("a24power")
    row <- g[g$family == "power_surface", ][1L, ]
    expect_equal(a24_arms_for(row), "lmm")
    expect_equal(length(a24_arms_for(g[g$family == "power_fig", ][1L, ])), 6L)
    g3 <- build_grid("a24sens")
    sa <- a24_arms_for(g3[1L, ])
    expect_equal(length(sa), 5L)
    expect_false("lmm_slope_only" %in% sa)
    row$rep_start <- 1L
    row$rep_end <- 2L
    res <- run_batch_slim(row, arms = a24_arms_for(row))
    expect_equal(nrow(res), 2L)
    expect_true(all(res$arm == "lmm"))
    row2 <- g[g$family == "power_fig", ][1L, ]
    row2$rep_start <- 1L
    row2$rep_end <- 1L
    expect_equal(nrow(run_batch_slim(row2, arms = a24_arms_for(row2))), 12L)
})
