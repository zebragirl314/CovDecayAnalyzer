## =========================================================================
## A2-4 (Tier 2) summaries. The per-replicate table is slim (see
## `run_batch_slim()`); cell parameters are joined from the grid with
## `a24_attach()`. Rejection is p_value <= alpha; a failed arm (NA p-value)
## counts as a non-rejection (Q6). Confidence intervals on single rates are
## Clopper-Pearson; paired differences use the Newcombe paired interval
## already in `paired_diff()`; unpaired differences use `.newcombe_diff()`
## with phi = 0.
## =========================================================================

SLIM_COLS <- c("scenario_id", "replicate_id", "bin_size", "arm",
    "test_variant", "status", "p_value", "F_stat", "df1", "df2", "singular",
    "conv_code", "sigma", "sd_slope", "sd_intercept", "re_corr",
    "p_equal_min", "slope_g1", "slope_g2", "slope_g3", "slope_g4",
    "n_samples", "n_bins", "n_warnings", "wall_time", "kw_perm_method",
    "dunn_perm_method")

#' Arms to run for one batch: the Stage 2 Sobol surface points need only the
#' primary LMM; Stage 3 sensitivity cells run every arm except the slope-only
#' comparator; every other cell runs all arms.
a24_arms_for <- function(batch) {
    if (all(batch$family == "power_surface")) return("lmm")
    ## Stage 3 drops the slope-only cautionary arm (see `a24_sens_grid()`).
    if (all(grepl("^sens_", batch$family))) {
        return(setdiff(names(ARM_FUNS), "lmm_slope_only"))
    }
    names(ARM_FUNS)
}

#' `run_batch()` with only the per-replicate result columns; scenario
#' parameters are recovered from the grid by `scenario_id`, which keeps the
#' Tier 2 result tables an order of magnitude smaller.
run_batch_slim <- function(batch, methods_files = NULL,
                           arms = names(ARM_FUNS)) {
    res <- run_batch(batch, methods_files, arms = arms)
    res[intersect(SLIM_COLS, names(res))]
}

## Vectorised Clopper-Pearson interval.
cp_vec <- function(x, n, level = 0.95) {
    a <- (1 - level) / 2
    lo <- ifelse(x == 0, 0, stats::qbeta(a, x, n - x + 1))
    hi <- ifelse(x == n, 1, stats::qbeta(1 - a, x + 1, n - x))
    list(lower = lo, upper = hi)
}

## Join cell parameters from the grid onto the slim replicate table.
CELL_COLS <- c("family", "sweep_factor", "sweep_level", "G", "n", "n_spec",
    "Delta_std", "cv", "s_ref", "sigma_atom", "L", "anchor_frac", "tau_eff",
    "delta_std", "seed_id", "sens_spec")
a24_attach <- function(sim, grid) {
    key <- paste(grid$scenario_id, grid$bin_size)
    stopifnot(!anyDuplicated(key))
    ix <- match(paste(sim$scenario_id, sim$bin_size), key)
    if (anyNA(ix)) stop("Replicate rows with no matching grid cell.")
    cols <- intersect(CELL_COLS, names(grid))
    cbind(sim, grid[ix, cols])
}

## Per (cell, arm, test_variant) counts; `extra_by` adds grouping columns
## carried from the cell.
a24_cell_counts <- function(sim, alpha = ALPHA, by_cell = "scenario_id") {
    dt <- data.table::as.data.table(sim)
    ok <- dt$status == "ok"
    ok[is.na(ok)] <- FALSE
    dt[, `:=`(.ok = ok)]
    out <- list()
    for (a in alpha) {
        dt[, .rej := !is.na(p_value) & p_value <= a]
        r <- dt[, list(n_datasets = .N, n_reject = sum(.rej),
            n_failed = sum(!.ok),
            singular_rate = if (any(.ok)) mean(as.numeric(singular[.ok]),
                na.rm = TRUE) else NA_real_),
            by = c(by_cell, "bin_size", "arm", "test_variant")]
        r$alpha <- a
        out[[length(out) + 1L]] <- r
    }
    out <- data.table::rbindlist(out)
    out$rate <- out$n_reject / out$n_datasets
    ci <- cp_vec(out$n_reject, out$n_datasets)
    out$ci_lower <- ci$lower
    out$ci_upper <- ci$upper
    out$failure_rate <- out$n_failed / out$n_datasets
    tibble::as_tibble(out)
}

## Add cell columns (once per scenario) to a counts table.
.add_cells <- function(counts, grid, cols = CELL_COLS) {
    ix <- match(paste(counts$scenario_id, counts$bin_size),
        paste(grid$scenario_id, grid$bin_size))
    cbind(counts, grid[ix, intersect(cols, names(grid))])
}

## Exact null rejection rate of chi-square KW at (G, n balanced or vector):
## enumeration when the assignment count is small, else Monte Carlo over
## 2e5 random permutations of the ranks.
kw_chisq_null_rate <- function(n, G, alpha = ALPHA, B = 2e5L, seed = 1L) {
    ng <- if (length(n) == 1L) rep(as.integer(n), G) else as.integer(n)
    N <- sum(ng)
    crit <- stats::qchisq(1 - alpha, length(ng) - 1L)
    n_assign <- exp(lgamma(N + 1) - sum(lgamma(ng + 1)))
    ranks <- seq_len(N)
    if (n_assign <= KW_PERM_MAX) {
        sums <- .enum_rank_sums(ranks, ng)
        H <- apply(sums, 1L, .kw_stat, ng = ng, N = N)
        mean(H >= crit - 1e-9)
    } else {
        g <- rep(seq_along(ng), ng)
        H <- with_seed(seed, vapply(seq_len(B), function(i) {
            .kw_stat(as.numeric(tapply(ranks, sample(g), sum)), ng, N)
        }, numeric(1L)))
        mean(H >= crit - 1e-9)
    }
}

.write_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(x, path, row.names = FALSE)
    path
}

## --------------------------- Stage 1: Type I ----------------------------

a24_typeI_main <- function(sim, grid) {
    cells <- grid[grid$family == "typeI_null", ]
    d <- sim[sim$scenario_id %in% cells$scenario_id, ]
    t <- .add_cells(a24_cell_counts(d, alpha = c(0.01, 0.05, 0.10)), grid)
    t$q7_class <- ifelse(t$alpha == 0.05,
        .q7_class(t$ci_lower, t$ci_upper), NA_character_)
    ## Exact null rate of chi-square KW beside the KW cells (A2-3 note).
    kw <- t$arm == "kw_chisq" & t$test_variant == "omnibus"
    key <- unique(t[kw, c("G", "n", "alpha")])
    ex <- vapply(seq_len(nrow(key)), function(i) {
        kw_chisq_null_rate(key$n[i], key$G[i], key$alpha[i])
    }, numeric(1L))
    t$kw_chisq_exact_null_rate <- NA_real_
    t$kw_chisq_exact_null_rate[kw] <- ex[match(
        paste(t$G, t$n, t$alpha)[kw], paste(key$G, key$n, key$alpha))]
    t <- t[order(t$arm, t$test_variant, t$G, t$n, t$alpha),
        c("arm", "test_variant", "G", "n", "alpha", "n_datasets", "n_reject",
            "rate", "ci_lower", "ci_upper", "q7_class", "failure_rate",
            "singular_rate", "kw_chisq_exact_null_rate", "scenario_id",
            "bin_size")]
    tibble::as_tibble(t)
}

## Family-wise error of the pairwise procedures among truly equal pairs:
## LMM Tukey (emmeans) ungated and gated on the omnibus p <= 0.05; Dunn-Holm
## (normal approximation and permutation) ungated and gated (the gated KW
## rows already carry the gate). Partial-null cells (one group differs) and
## the G >= 3 global-null cells.
a24_typeI_posthoc <- function(sim, grid) {
    cells <- grid[grid$family %in% c("typeI_partial", "typeI_null") &
        grid$G >= 3L, ]
    d <- data.table::as.data.table(sim[sim$scenario_id %in% cells$scenario_id, ])
    omni <- d[d$arm == "lmm" & d$test_variant == "omnibus",
        list(scenario_id, bin_size, replicate_id, status, omni_p = p_value,
            p_equal_min)]
    lmm_u <- data.table::copy(omni)
    lmm_u[, `:=`(method = "lmm_tukey_ungated",
        rej = !is.na(p_equal_min) & p_equal_min <= ALPHA)]
    lmm_g <- data.table::copy(omni)
    lmm_g[, `:=`(method = "lmm_tukey_gated",
        rej = !is.na(p_equal_min) & p_equal_min <= ALPHA &
            !is.na(omni_p) & omni_p <= ALPHA)]
    dunn <- d[d$test_variant %in% c("dunn_ungated", "dunn_gated",
        "dunn_perm_ungated", "dunn_perm_gated"),
        list(scenario_id, bin_size, replicate_id, status,
            method = paste(arm, test_variant, sep = "_"),
            rej = !is.na(p_equal_min) & p_equal_min <= ALPHA)]
    all <- data.table::rbindlist(list(
        lmm_u[, list(scenario_id, bin_size, replicate_id, status, method, rej)],
        lmm_g[, list(scenario_id, bin_size, replicate_id, status, method, rej)],
        dunn))
    r <- all[, list(n_datasets = .N, n_false_reject = sum(rej),
        failure_rate = mean(status != "ok")),
        by = list(scenario_id, bin_size, method)]
    r$fwer <- r$n_false_reject / r$n_datasets
    ci <- cp_vec(r$n_false_reject, r$n_datasets)
    r$ci_lower <- ci$lower
    r$ci_upper <- ci$upper
    r <- .add_cells(tibble::as_tibble(r), grid)
    r$cell_kind <- ifelse(r$family == "typeI_partial", "partial_null",
        "global_null")
    r <- r[order(r$method, r$cell_kind, r$G, r$n),
        c("method", "cell_kind", "G", "n", "Delta_std", "n_datasets",
            "n_false_reject", "fwer", "ci_lower", "ci_upper",
            "failure_rate", "scenario_id", "bin_size")]
    tibble::as_tibble(r)
}

write_a24_typeI <- function(sim, grid, out_dir = "results/tier2/typeI") {
    dir.create(file.path(out_dir, "raw"), recursive = TRUE,
        showWarnings = FALSE)
    paths <- c(
        .write_csv(a24_typeI_main(sim, grid),
            file.path(out_dir, "tbl-typeI-main.csv")),
        .write_csv(a24_typeI_posthoc(sim, grid),
            file.path(out_dir, "tbl-typeI-posthoc.csv")),
        .write_csv(as.data.frame(grid[setdiff(names(grid), "reps")]),
            file.path(out_dir, "grid.csv")))
    raw <- file.path(out_dir, "raw", "tier2-typeI.rds")
    saveRDS(sim, raw, compress = "xz")
    c(paths, raw)
}

## ---------------------------- Stage 2: power ----------------------------

a24_power_main <- function(sim, grid) {
    cells <- grid[grid$family == "power_fig", ]
    d <- sim[sim$scenario_id %in% cells$scenario_id, ]
    t <- .add_cells(a24_cell_counts(d, alpha = ALPHA), grid)
    t <- t[order(t$arm, t$test_variant, t$G, t$Delta_std, t$n),
        c("arm", "test_variant", "G", "n", "Delta_std", "cv", "n_datasets",
            "n_reject", "rate", "ci_lower", "ci_upper", "failure_rate",
            "singular_rate", "scenario_id", "bin_size")]
    tibble::as_tibble(t)
}

## Smallest n on the figure grid with power >= 0.8, per (arm omnibus, G,
## Delta_std, cv); also the smallest n whose CP lower bound is >= 0.8.
## The figure grid has one cv (the reference), so cv is constant.
a24_n_for_80 <- function(power_main) {
    d <- power_main[power_main$test_variant == "omnibus", ]
    key <- unique(d[c("arm", "G", "Delta_std", "cv")])
    out <- lapply(seq_len(nrow(key)), function(i) {
        k <- key[i, ]
        x <- d[d$arm == k$arm & d$G == k$G & d$Delta_std == k$Delta_std &
            d$cv == k$cv, ]
        x <- x[order(x$n), ]
        first <- function(ix) if (any(ix)) x$n[which(ix)[1L]] else NA_integer_
        cbind(k, tibble::tibble(
            n_for_80pct = first(x$rate >= 0.8),
            n_for_80pct_lower_ci = first(x$ci_lower >= 0.8),
            max_n_tested = max(x$n),
            power_at_max_n = x$rate[which.max(x$n)]))
    })
    t <- do.call(rbind, out)
    t[order(t$arm, t$G, t$Delta_std), ]
}

a24_surface_data <- function(sim, grid) {
    cells <- grid[grid$family == "power_surface", ]
    d <- sim[sim$scenario_id %in% cells$scenario_id &
        sim$test_variant == "omnibus", ]
    t <- .add_cells(a24_cell_counts(d, alpha = ALPHA), grid)
    t$point <- as.integer(t$sweep_level)
    t <- t[order(t$point, t$arm),
        c("point", "arm", "n", "G", "Delta_std", "cv", "sigma_atom", "L",
            "anchor_frac", "tau_eff", "delta_std", "n_datasets", "n_reject",
            "rate", "ci_lower", "ci_upper", "failure_rate", "singular_rate",
            "scenario_id")]
    tibble::as_tibble(t)
}

## Collapse check (PLAN.md 2.5): power at matched Delta_std / cv (cv fixed
## at 1/3 here) across s_ref in {3e-6, 6e-6, 1.2e-5}, G = 2, n in {3, 5,
## 10}, Delta_std in {0.5, 1}. The s_ref = 6e-6 cells are the figure-grid
## cells. Datasets are independent across s_ref (different scenarios), so
## the difference interval is unpaired Newcombe.
a24_collapse_table <- function(sim, grid, arms = c("lmm", "two_stage_pooled")) {
    cells <- grid[(grid$family == "power_collapse") |
        (grid$family == "power_fig" & grid$G == 2L &
            grid$n %in% c(3L, 5L, 10L) & grid$Delta_std %in% c(0.5, 1)), ]
    d <- sim[sim$scenario_id %in% cells$scenario_id &
        sim$test_variant == "omnibus" & sim$arm %in% arms, ]
    t <- .add_cells(a24_cell_counts(d, alpha = ALPHA), grid)
    out <- list()
    for (arm in arms) for (n in c(3L, 5L, 10L)) for (dl in c(0.5, 1)) {
        x <- t[t$arm == arm & t$n == n & t$Delta_std == dl, ]
        x <- x[order(x$s_ref), ]
        if (nrow(x) < 2L) next
        for (i in seq_len(nrow(x) - 1L)) for (j in (i + 1L):nrow(x)) {
            nd <- .newcombe_diff(x$n_reject[j], x$n_datasets[j],
                x$n_reject[i], x$n_datasets[i])
            out[[length(out) + 1L]] <- tibble::tibble(arm = arm, n = n,
                Delta_std = dl, cv = x$cv[i], s_ref_a = x$s_ref[i],
                s_ref_b = x$s_ref[j], rate_a = x$rate[i], rate_b = x$rate[j],
                n_a = x$n_datasets[i], n_b = x$n_datasets[j],
                diff_b_minus_a = nd[["diff"]], ci_lower = nd[["lower"]],
                ci_upper = nd[["upper"]],
                exceeds_0.03 = abs(nd[["diff"]]) > 0.03,
                ci_excludes_0 = nd[["lower"]] > 0 | nd[["upper"]] < 0)
        }
    }
    do.call(rbind, out)
}

write_a24_power <- function(sim, grid, out_dir = "results/tier2/power") {
    dir.create(file.path(out_dir, "raw"), recursive = TRUE,
        showWarnings = FALSE)
    pm <- a24_power_main(sim, grid)
    paths <- c(
        .write_csv(pm, file.path(out_dir, "tbl-power-main.csv")),
        .write_csv(a24_n_for_80(pm), file.path(out_dir, "tbl-n-for-80pct.csv")),
        .write_csv(a24_surface_data(sim, grid),
            file.path(out_dir, "surface", "power-surface-data.csv")),
        .write_csv(a24_collapse_table(sim, grid),
            file.path(out_dir, "tbl-collapse-tier2.csv")),
        .write_csv(as.data.frame(grid[setdiff(names(grid), "reps")]),
            file.path(out_dir, "grid.csv")))
    raw <- file.path(out_dir, "raw", "tier2-power.rds")
    saveRDS(sim, raw, compress = "xz")
    c(paths, raw)
}

## --------------------------- Stage 3: sensitivities ----------------------

SENS_FILES <- c(
    sens_S0 = "tbl-sens-unbalanced.csv",
    sens_S1a = "tbl-sens-autocorr.csv",
    sens_S1b = "tbl-sens-misspec.csv",
    sens_S2 = "tbl-sens-dropout.csv",
    sens_S3 = "tbl-sens-hetero-residual.csv",
    sens_S4 = "tbl-sens-hetero-slope.csv",
    sens_S5 = "tbl-sens-anchor-extreme.csv")
SENS_PAIRED <- c("sens_S1a", "sens_S1b", "sens_S2", "sens_S3", "sens_S4")

## Rejection indicators (failure = non-rejection) by replicate for one
## (scenario, arm, variant): a numeric/logical vector ordered by replicate.
.rej_index <- function(sim) {
    dt <- data.table::as.data.table(sim)
    dt[, .rej := !is.na(p_value) & p_value <= ALPHA]
    data.table::setkey(dt, scenario_id, arm, test_variant, replicate_id)
    dt
}

.get_rej <- function(dt, qid, qarm, qvar) {
    x <- dt[list(qid, qarm, qvar)]
    data.table::setorder(x, replicate_id)
    list(rej = x$.rej, rep = x$replicate_id, ok = x$status == "ok")
}

## Reference scenario for each sensitivity cell. Paired families share the
## baseline draws (seed_id = reference scenario_id); S5 compares with the
## reference cell of the same G, n, Delta_std (different geometry, unpaired);
## S0 compares with the balanced cell whose n equals the smaller group's n
## when that is a reference n (3 or 5), else none.
.sens_reference <- function(cell, grid) {
    ref <- grid[grid$family == "sens_ref", ]
    if (cell$family %in% SENS_PAIRED) {
        return(list(id = cell$seed_id, paired = TRUE))
    }
    pick <- function(n) {
        r <- ref[ref$G == cell$G & ref$n == n & ref$Delta_std == cell$Delta_std, ]
        if (nrow(r) == 1L) r$scenario_id else NA_character_
    }
    if (cell$family == "sens_S5") return(list(id = pick(cell$n), paired = FALSE))
    if (cell$family == "sens_S0") {
        nn <- min(as.integer(strsplit(cell$n_spec, ",")[[1L]]))
        return(list(id = if (nn %in% c(3L, 5L)) pick(nn) else NA_character_,
            paired = FALSE))
    }
    list(id = NA_character_, paired = FALSE)
}

## One sensitivity table. Per (cell, arm, variant): the sensitivity rate with
## CP CI; the reference rate with CP CI; their difference with a CI (paired
## Newcombe when the datasets share seeds, else unpaired); metric is "typeI"
## when Delta_std = 0 and "power" otherwise.
a24_sens_table <- function(sim, grid, family) {
    cells <- grid[grid$family == family, ]
    dt <- .rej_index(sim)
    av <- unique(dt[, list(arm, test_variant)])
    rows <- list()
    for (i in seq_len(nrow(cells))) {
        cell <- cells[i, ]
        ref <- .sens_reference(cell, grid)
        for (k in seq_len(nrow(av))) {
            a <- av$arm[k]
            v <- av$test_variant[k]
            s <- .get_rej(dt, cell$scenario_id, a, v)
            x <- sum(s$rej)
            nn <- length(s$rej)
            ci <- cp_vec(x, nn)
            r <- tibble::tibble(
                family = family, sweep_factor = cell$sweep_factor,
                sweep_level = cell$sweep_level, G = cell$G,
                n_spec = cell$n_spec, Delta_std = cell$Delta_std,
                metric = if (cell$Delta_std == 0) "typeI" else "power",
                arm = a, test_variant = v, n_datasets = nn,
                rate = x / nn, ci_lower = ci$lower, ci_upper = ci$upper,
                failure_rate = mean(!s$ok),
                ref_scenario_id = ref$id, pairing = NA_character_,
                ref_rate = NA_real_, ref_ci_lower = NA_real_,
                ref_ci_upper = NA_real_, diff = NA_real_,
                diff_ci_lower = NA_real_, diff_ci_upper = NA_real_,
                sens_only_rejects = NA_real_, ref_only_rejects = NA_real_)
            if (!is.na(ref$id)) {
                b <- .get_rej(dt, ref$id, a, v)
                xr <- sum(b$rej)
                nr <- length(b$rej)
                cir <- cp_vec(xr, nr)
                r$ref_rate <- xr / nr
                r$ref_ci_lower <- cir$lower
                r$ref_ci_upper <- cir$upper
                if (ref$paired) {
                    m <- intersect(s$rep, b$rep)
                    pd <- paired_diff(s$rej[match(m, s$rep)],
                        b$rej[match(m, b$rep)])
                    r$pairing <- "paired"
                    r$diff <- pd$diff
                    r$diff_ci_lower <- pd$ci_lower
                    r$diff_ci_upper <- pd$ci_upper
                    r$sens_only_rejects <- pd$a_only
                    r$ref_only_rejects <- pd$b_only
                } else {
                    nd <- .newcombe_diff(x, nn, xr, nr)
                    r$pairing <- "unpaired"
                    r$diff <- nd[["diff"]]
                    r$diff_ci_lower <- nd[["lower"]]
                    r$diff_ci_upper <- nd[["upper"]]
                }
            }
            rows[[length(rows) + 1L]] <- r
        }
    }
    do.call(rbind, rows)
}

## Paired primary-LMM minus pooled-two-stage rejection rate at every cell of
## the sensitivity stage (the paper's methodological comparison). Positive =
## LMM rejects more. The same datasets feed both arms.
a24_lmm_vs_pooled <- function(sim, grid, comparator = "two_stage_pooled") {
    dt <- .rej_index(sim)
    cells <- grid[grepl("^sens_", grid$family), ]
    rows <- lapply(seq_len(nrow(cells)), function(i) {
        cell <- cells[i, ]
        a <- .get_rej(dt, cell$scenario_id, "lmm", "omnibus")
        b <- .get_rej(dt, cell$scenario_id, comparator, "omnibus")
        pd <- paired_diff(a$rej, b$rej)
        cbind(tibble::tibble(family = cell$family,
            sweep_factor = cell$sweep_factor, sweep_level = cell$sweep_level,
            G = cell$G, n_spec = cell$n_spec, Delta_std = cell$Delta_std,
            metric = if (cell$Delta_std == 0) "typeI" else "power",
            comparator = comparator), pd)
    })
    t <- do.call(rbind, rows)
    names(t)[names(t) == "rate_a"] <- "rate_lmm"
    names(t)[names(t) == "rate_b"] <- "rate_comparator"
    names(t)[names(t) == "a_only"] <- "lmm_only_rejects"
    names(t)[names(t) == "b_only"] <- "comparator_only_rejects"
    t
}

## S3 rerun table: one row per cell with the LMM, pooled two-stage and Welch
## two-stage omnibus rates (CP CIs), the paired LMM-minus-pooled difference
## (same datasets), and tau_eff / tau for the cell.
write_a24_s3rerun <- function(sim, grid, out_dir = "results/tier2/sensitivities",
                              stem = "hetero-residual-rerun", grid_file = "grid-s3rerun.csv",
                              raw_file = "tier2-s3rerun.rds") {
    dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
    dt <- .rej_index(sim)
    rows <- lapply(seq_len(nrow(grid)), function(i) {
        cell <- grid[i, ]
        g <- function(a) .get_rej(dt, cell$scenario_id, a, "omnibus")
        l <- g("lmm"); p <- g("two_stage_pooled"); w <- g("two_stage")
        one <- function(z, nm) {
            ci <- cp_vec(sum(z$rej), length(z$rej))
            out <- tibble::tibble(rate = sum(z$rej) / length(z$rej),
                ci_lower = ci$lower, ci_upper = ci$upper,
                failure_rate = mean(!z$ok))
            names(out) <- paste0(names(out), "_", nm)
            out
        }
        pd <- paired_diff(l$rej, p$rej)
        pw <- paired_diff(l$rej, w$rej)
        cbind(tibble::tibble(family = cell$family, sigma_atom = cell$sigma_atom,
            G = cell$G, n = cell$n, Delta_std = cell$Delta_std,
            metric = if (cell$Delta_std == 0) "typeI" else "power",
            n_datasets = length(l$rej), tau_eff_over_tau = cell$tau_eff / cell$tau),
            one(l, "lmm"), one(p, "pooled"), one(w, "welch"),
            tibble::tibble(lmm_minus_pooled = pd$diff,
                lmm_minus_pooled_ci_lower = pd$ci_lower,
                lmm_minus_pooled_ci_upper = pd$ci_upper,
                lmm_only_rejects = pd$a_only, pooled_only_rejects = pd$b_only,
                mcnemar_p = pd$mcnemar_p,
                lmm_minus_welch = pw$diff,
                lmm_minus_welch_ci_lower = pw$ci_lower,
                lmm_minus_welch_ci_upper = pw$ci_upper))
    })
    t <- do.call(rbind, rows)
    paths <- c(
        .write_csv(t, file.path(out_dir, paste0("tbl-sens-", stem, ".csv"))),
        .write_csv(as.data.frame(grid[setdiff(names(grid), "reps")]),
            file.path(out_dir, grid_file)))
    raw <- file.path(out_dir, "raw", raw_file)
    saveRDS(sim, raw, compress = "xz")
    unname(c(paths, raw))
}

write_a24_s3baseline <- function(sim, grid, out_dir = "results/tier2/sensitivities") {
    write_a24_s3rerun(sim, grid, out_dir, stem = "hetero-residual-baseline",
        grid_file = "grid-s3baseline.csv", raw_file = "tier2-s3baseline.rds")
}

write_a24_sens <- function(sim, grid, out_dir = "results/tier2/sensitivities") {
    dir.create(file.path(out_dir, "raw"), recursive = TRUE,
        showWarnings = FALSE)
    paths <- vapply(names(SENS_FILES), function(f) {
        .write_csv(a24_sens_table(sim, grid, f),
            file.path(out_dir, SENS_FILES[[f]]))
    }, character(1L))
    paths <- c(paths,
        .write_csv(a24_lmm_vs_pooled(sim, grid),
            file.path(out_dir, "tbl-sens-lmm-vs-pooled.csv")),
        .write_csv(as.data.frame(grid[setdiff(names(grid), "reps")]),
            file.path(out_dir, "grid.csv")))
    raw <- file.path(out_dir, "raw", "tier2-sensitivities.rds")
    saveRDS(sim, raw, compress = "xz")
    unname(c(paths, raw))
}
