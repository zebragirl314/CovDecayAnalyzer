## =========================================================================
## Generator (PLAN.md 2.1 and 2.2): simulate BEFORE normalization.
##
##   1. Draw per-sample mean coverage on a 1 kb atom grid.
##   2. Aggregate atoms to the analysis bin size with the same arithmetic as
##      the package's `bin_coverage()`.
##   3. Call the package's own `normalize_coverage()` per sample.
##
## Baseline model, sample i in group g(i), atom a at midpoint x_a:
##   mean_cov_ia = D_i * exp(-s_i * |x_a - anchor| + e_ia)
##   s_i  = s_ref + Delta_g(i) + u_i,  u_i ~ N(0, tau^2)
##   D_i  ~ LogNormal(log depth_median, depth_sd_log^2)
##   e_ia ~ N(0, sigma_atom^2), independent
## Group pattern (Q9): one group (the last) has s_ref + Delta; the rest
## share s_ref.
## =========================================================================

## 1 kb atom grid for a genome of length L. The last atom is truncated at L.
atom_geometry <- function(L, atom_size = ATOM_SIZE) {
    n_atoms <- as.integer(ceiling(L / atom_size))
    start <- (seq_len(n_atoms) - 1L) * atom_size + 1L
    end <- pmin(start + atom_size - 1L, L)
    data.frame(atom_start = start, atom_end = end,
        atom_mid = (start + end) / 2, width = end - start + 1L)
}

## Bin geometry using the same formulas as `bin_coverage()`:
## bin_start = (bin - 1) * bin_size + 1; bin_end = pmin(start + size - 1, L);
## bin_mid = (start + end) / 2.
bin_geometry <- function(L, bin_size, atom_size = ATOM_SIZE) {
    n_bins <- as.integer(ceiling(L / bin_size))
    bin_start <- (seq_len(n_bins) - 1L) * bin_size + 1L
    bin_end <- pmin(bin_start + bin_size - 1L, L)
    data.frame(bin = seq_len(n_bins), bin_start = bin_start,
        bin_end = bin_end, bin_mid = (bin_start + bin_end) / 2)
}

## Per-sample design and truth: group, slope, depth. RNG order is fixed
## (u, then D, then atom noise) so every component is stable across bin
## sizes and across sensitivities that add draws afterwards.
.draw_samples <- function(sc) {
    npg <- n_per_group(sc)
    group <- factor(rep(paste0("g", seq_len(sc$G)), times = npg),
        levels = paste0("g", seq_len(sc$G)))
    N <- length(group)
    delta <- group_deltas(sc$G, sc$Delta)[as.integer(group)]
    sample <- paste0(group, "_s", ave(seq_len(N), group, FUN = seq_along))
    samples <- tibble::tibble(
        sample = sample, group = group, delta = delta,
        tau_i = rep(sc$tau, N), sigma_i = rep(sc$sigma_atom, N))
    if (sens_active(sc, "S3")) {
        samples$s_group <- sc$s_ref + delta
        samples <- sens_s3_scales(samples, sc)
    }
    u <- stats::rnorm(N) * samples$tau_i
    D <- stats::rlnorm(N, meanlog = log(sc$depth_median),
        sdlog = sc$depth_sd_log)
    samples$u_i <- u
    samples$s_i <- sc$s_ref + delta + u
    samples$D <- D
    samples
}

## Atom-level mean coverage (atoms x samples) with attributes for truth.
simulate_atoms <- function(sc, overrides = NULL) {
    samples <- .draw_samples(sc)
    if (!is.null(overrides$D)) {
        if (length(overrides$D) != nrow(samples)) {
            stop("`overrides$D` must have one value per sample.")
        }
        samples$D <- overrides$D
    }
    atoms <- atom_geometry(sc$L, sc$atom_size)
    N <- nrow(samples)
    n_atoms <- nrow(atoms)
    dist <- abs(atoms$atom_mid - sc$anchor)

    ## Base noise is always drawn, so streams are shared with sensitivities.
    e <- matrix(stats::rnorm(n_atoms * N), n_atoms, N)
    e <- sweep(e, 2L, samples$sigma_i, `*`)

    if (sens_active(sc, "S1a")) {
        e <- sens_s1a_noise(e, sc)
    }
    log_mean <- if (sens_active(sc, "S1b")) {
        sens_s1b_curve(dist, samples$s_i, sc)
    } else {
        -outer(dist, samples$s_i)
    }
    cov <- sweep(exp(log_mean + e), 2L, samples$D, `*`)

    if (sens_active(sc, "S4")) cov <- sens_s4_counts(cov, sc)
    if (sens_active(sc, "S2")) cov <- sens_s2_dropout(cov, sc)

    list(cov = cov, atoms = atoms, samples = samples)
}

## Aggregate atoms to bins: width-weighted mean over the atoms in a bin,
## which equals `bin_coverage()`'s per-base mean when coverage is constant
## within an atom (it differs from an unweighted atom mean only for a
## truncated final atom).
bin_atoms <- function(atom_cov, atoms, L, bin_size) {
    if (bin_size %% ATOM_SIZE != 0L) {
        stop("`bin_size` must be a multiple of the atom size.")
    }
    bin <- (atoms$atom_start - 1L) %/% bin_size + 1L
    num <- rowsum(atom_cov * atoms$width, group = bin, reorder = TRUE)
    den <- rowsum(atoms$width, group = bin, reorder = TRUE)
    mean_cov <- num / as.vector(den)
    geo <- bin_geometry(L, bin_size)
    list(mean_cov = mean_cov, geo = geo)
}

#' Generate one replicate dataset.
#'
#' @param scenario A `cda_scenario` from [make_scenario()].
#' @param replicate_id Integer replicate index; with the scenario ID it
#'   fixes the seed.
#' @param bin_size Analysis bin size in bp; defaults to the scenario's.
#' @param overrides Optional list; `D` replaces the drawn depths (used by
#'   the normalization check).
#' @return A tibble with columns `sample`, `group`, `bin`, `bin_start`,
#'   `bin_end`, `bin_mid`, `mean_cov`, `norm_cov`, `log_norm_cov`,
#'   `distance_from_anchor`, ready for `test_decay_slopes()`. Attributes
#'   `scenario_id`, `replicate_id`, `seed`, `bin_size`, and `truth` (a
#'   tibble of per-sample generating values) are attached.
generate_replicate <- function(scenario, replicate_id,
                               bin_size = scenario$bin_size,
                               overrides = NULL) {
    ensure_package()
    seed <- derive_seed(scenario$scenario_id, replicate_id)
    sim <- with_seed(seed, simulate_atoms(scenario, overrides))

    binned <- bin_atoms(sim$cov, sim$atoms, scenario$L, bin_size)
    geo <- binned$geo
    normalize <- get("normalize_coverage", envir = asNamespace("CovDecayAnalyzer"))

    per_sample <- lapply(seq_len(nrow(sim$samples)), function(i) {
        dt <- data.table::data.table(
            bin = geo$bin, mean_cov = binned$mean_cov[, i],
            bin_start = geo$bin_start, bin_end = geo$bin_end,
            bin_mid = geo$bin_mid)
        dt <- normalize(dt)
        dt[, sample := sim$samples$sample[i]]
        dt[, group := as.character(sim$samples$group[i])]
        dt
    })
    out <- data.table::rbindlist(per_sample)

    anchor <- if (sens_active(scenario, "S5")) {
        sens_s5_analysis_anchor(scenario)
    } else {
        scenario$anchor
    }
    out[, distance_from_anchor := abs(bin_mid - anchor)]
    data.table::setcolorder(out, c("sample", "group", "bin", "bin_start",
        "bin_end", "bin_mid", "mean_cov", "norm_cov", "log_norm_cov",
        "distance_from_anchor"))
    out <- tibble::as_tibble(out)
    out$group <- factor(out$group, levels = levels(sim$samples$group))

    attr(out, "scenario_id") <- scenario$scenario_id
    attr(out, "replicate_id") <- replicate_id
    attr(out, "seed") <- seed
    attr(out, "bin_size") <- bin_size
    attr(out, "analysis_anchor") <- anchor
    attr(out, "truth") <- sim$samples
    out
}
