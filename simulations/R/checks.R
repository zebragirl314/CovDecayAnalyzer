## =========================================================================
## Generator validation checks 1 to 6 (PLAN.md section 6).
##
## Each check returns list(name, pass, detail). The testthat file wraps them
## in `test_that()` blocks and `scripts/check_generator.R` prints them, so
## the two cannot drift apart. These check the instrument (the generator
## and the analysis wrapper), not the claims the study makes about the test.
##
## Replicate counts default to values that finish in minutes and can be
## raised with CDA_SIM_CHECK_R (calibration) for a tighter check.
## =========================================================================

.check <- function(name, pass, detail) {
    list(name = name, pass = isTRUE(pass), detail = detail)
}

.check_R <- function(default) {
    as.integer(Sys.getenv("CDA_SIM_CHECK_R", as.character(default)))
}

## 1. Parameter recovery: at large n the fitted Delta, tau, and sigma
## recover the generating values.
check_parameter_recovery <- function(n = 40L, reps = 8L, Delta = 3e-6) {
    sc <- make_scenario(G = 2L, n = n, Delta = Delta)
    rows <- do.call(rbind, lapply(seq_len(reps), function(r) {
        analyze_lmm(generate_replicate(sc, r))
    }))
    ok <- all(rows$status == "ok")
    delta_hat <- rows$slope_g1 - rows$slope_g2  # LMM slope is -s
    se_mean <- sc$tau * sqrt(2 / n) / sqrt(reps)
    sigma_expected <- sc$sigma_atom / sqrt(sc$bin_size / sc$atom_size)
    d_err <- mean(delta_hat) - Delta
    tau_rel <- mean(rows$sd_slope) / sc$tau - 1
    sig_rel <- mean(rows$sigma) / sigma_expected - 1
    pass <- ok && abs(d_err) < 4 * se_mean && abs(tau_rel) < 0.12 &&
        abs(sig_rel) < 0.05
    .check("parameter recovery", pass, sprintf(
        paste("n=%d/group, %d reps, all fits ok=%s. Delta: mean %.3g vs",
            "%.3g (error %.2g, 4 SE = %.2g). tau: mean %.3g vs %.3g",
            "(%+.1f%%). sigma: mean %.4f vs %.4f (%+.1f%%)."),
        n, reps, ok, mean(delta_hat), Delta, d_err, 4 * se_mean,
        mean(rows$sd_slope), sc$tau, 100 * tau_rel,
        mean(rows$sigma), sigma_expected, 100 * sig_rel))
}

## 2. Idealized calibration: with the correctly specified baseline model
## and large n, Type I is near 0.05 and the p-value ECDF is near uniform.
check_calibration <- function(n = 15L, R = .check_R(200L)) {
    sc <- make_scenario(G = 2L, n = n, Delta = 0)
    rows <- do.call(rbind, lapply(seq_len(R), function(r) {
        analyze_lmm(generate_replicate(sc, r))
    }))
    p <- rows$p_value[rows$status == "ok"]
    rej <- sum(p <= 0.05)
    ci <- cp_interval(rej, length(p))
    ks <- suppressWarnings(stats::ks.test(p, "punif"))
    pass <- length(p) == R && ci[["lower"]] <= 0.05 && ci[["upper"]] >= 0.05 &&
        ks$p.value > 0.01
    .check("idealized calibration", pass, sprintf(
        paste("n=%d/group, R=%d, fits ok=%d. Type I at 0.05 = %.3f",
            "(Clopper-Pearson 95%% CI %.3f to %.3f). KS test of uniformity",
            "p=%.3f."),
        n, R, length(p), rej / length(p), ci[["lower"]], ci[["upper"]],
        ks$p.value))
}

## 3. Positive control: a very large Delta gives power near 1.
check_positive_control <- function(n = 4L, R = .check_R(30L),
                                   Delta = 3e-5) {
    sc <- make_scenario(G = 2L, n = n, Delta = Delta)
    rows <- do.call(rbind, lapply(seq_len(R), function(r) {
        analyze_lmm(generate_replicate(sc, r))
    }))
    p <- rows$p_value[rows$status == "ok"]
    power <- mean(p <= 0.05)
    pass <- length(p) == R && power >= 0.95
    .check("positive control", pass, sprintf(
        "n=%d/group, Delta=%.1e (Delta/tau=%.0f), R=%d: power %.3f.",
        n, Delta, Delta / sc$tau, R, power))
}

## 4. Path equivalence: the atom path and the package's per-base
## `load_all_samples()` path agree on `log_norm_cov`.
check_path_equivalence <- function(bin_sizes = c(5000L, 1000L),
                                   L = 238080) {
    ensure_package()
    sc <- make_scenario(G = 2L, n = 2L, Delta = 3e-6, L = L)
    seed <- derive_seed(sc$scenario_id, 1L)
    sim <- with_seed(seed, simulate_atoms(sc))
    dir <- tempfile("cda-perbase-")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)

    files <- vapply(seq_len(nrow(sim$samples)), function(i) {
        f <- file.path(dir, paste0(sim$samples$sample[i], ".tab"))
        cov <- rep(sim$cov[, i], times = sim$atoms$width)
        data.table::fwrite(data.table::data.table(
            position = seq_len(L), unique_top_cov = cov,
            unique_bot_cov = 0, redundant_top_cov = 0,
            redundant_bot_cov = 0), f)  # 15 significant digits by default
        f
    }, character(1L))
    meta <- data.frame(sample = sim$samples$sample,
        group = as.character(sim$samples$group), coverage_file = files)

    diffs <- vapply(bin_sizes, function(b) {
        atom_path <- generate_replicate(sc, 1L, bin_size = b)
        base_path <- CovDecayAnalyzer::load_all_samples(meta, bin_size = b,
            genome_length = as.integer(L))
        a <- atom_path[order(atom_path$sample, atom_path$bin), ]
        p <- as.data.frame(base_path)
        p <- p[order(p$sample, p$bin), ]
        stopifnot(identical(a$sample, p$sample), identical(a$bin, p$bin))
        max(abs(a$log_norm_cov - p$log_norm_cov),
            abs(a$bin_mid - p$bin_mid))
    }, numeric(1L))
    pass <- all(diffs < 1e-8)
    .check("path equivalence", pass, sprintf(
        "L=%g (final atom truncated: %s). Max abs difference in log_norm_cov and bin_mid: %s.",
        L, L %% ATOM_SIZE != 0,
        paste(sprintf("%.2e at %d bp", diffs, bin_sizes), collapse = "; ")))
}

## 5. Determinism: a replicate regenerates exactly from its seed, adding or
## reordering replicates leaves existing ones unchanged, and generating
## does not disturb the caller's RNG stream.
check_determinism <- function() {
    sc <- make_scenario(G = 3L, n = 3L, Delta = 2e-6)
    alone <- generate_replicate(sc, 3L)
    invisible(generate_replicate(sc, 1L))
    invisible(generate_replicate(sc, 2L))
    after_others <- generate_replicate(sc, 3L)
    same_data <- identical(as.data.frame(alone), as.data.frame(after_others))

    r1 <- analyze_lmm(alone)
    r2 <- analyze_lmm(after_others)
    same_p <- identical(r1$p_value, r2$p_value)

    set.seed(99L)
    expected <- stats::runif(1L)
    set.seed(99L)
    invisible(generate_replicate(sc, 4L))
    rng_untouched <- identical(stats::runif(1L), expected)

    other <- generate_replicate(sc, 4L)
    differs <- !identical(as.data.frame(alone), as.data.frame(other))

    ## Paired bin sizes share the same simulated samples.
    t5 <- attr(generate_replicate(sc, 5L, bin_size = 5000L), "truth")
    t1 <- attr(generate_replicate(sc, 5L, bin_size = 1000L), "truth")
    paired <- identical(t5, t1)

    pass <- same_data && same_p && rng_untouched && differs && paired
    .check("determinism", pass, sprintf(
        paste("regenerated replicate identical: %s; p-value identical: %s;",
            "caller RNG untouched: %s; different replicates differ: %s;",
            "5 kb and 1 kb share truth: %s."),
        same_data, same_p, rng_untouched, differs, paired))
}

## 6. Normalization is exercised: changing depth D_i alone changes the raw
## coverage but no normalized value and no test result.
check_normalization <- function() {
    sc <- make_scenario(G = 2L, n = 4L, Delta = 3e-6)
    base <- generate_replicate(sc, 1L)
    D <- attr(base, "truth")$D
    scaled <- generate_replicate(sc, 1L,
        overrides = list(D = D * c(0.2, 5, 1.7, 0.4, 3, 0.9, 2.2, 0.6)))
    raw_changed <- max(abs(base$mean_cov - scaled$mean_cov)) > 1
    norm_diff <- max(abs(base$log_norm_cov - scaled$log_norm_cov))
    pb <- analyze_lmm(base)$p_value
    ps <- analyze_lmm(scaled)$p_value
    p_diff <- abs(pb - ps)
    ## The normalized data agree to rounding error. The fitted p-value agrees
    ## only to optimizer tolerance (observed about 5e-6): the fits sit on the
    ## singular boundary (M1), where REML is flat and the optimizer stops at
    ## slightly different points for data that differ at 1e-16.
    pass <- raw_changed && norm_diff < 1e-9 && p_diff < 1e-4
    .check("normalization exercised", pass, sprintf(
        paste("raw mean_cov changed: %s; max |change| in log_norm_cov:",
            "%.2e; |change| in p-value: %.2e."),
        raw_changed, norm_diff, p_diff))
}

run_all_checks <- function() {
    list(
        check_parameter_recovery(),
        check_calibration(),
        check_positive_control(),
        check_path_equivalence(),
        check_determinism(),
        check_normalization()
    )
}
