## =========================================================================
## Summaries. A2-1 defines only what Tier 0 needs; the Type I and power
## summaries (Clopper-Pearson, Q7) are completed in A2-2.
## =========================================================================

## Exact (Clopper-Pearson) interval for a proportion.
cp_interval <- function(x, n, level = 0.95) {
    if (n == 0L) return(c(lower = NA_real_, upper = NA_real_))
    ci <- stats::binom.test(x, n, conf.level = level)$conf.int
    c(lower = ci[1L], upper = ci[2L])
}

#' Wall-time summary by scenario and bin size.
#'
#' @param sim Tibble of replicate rows from [run_batch()].
#' @return Tibble: one row per scenario and bin size with the number of
#'   fits, the number of bins and samples, and the median, 95th percentile,
#'   and maximum wall time in seconds.
summarize_timing <- function(sim) {
    if ("arm" %in% names(sim)) sim <- sim[sim$arm == "lmm", ]
    key <- paste(sim$scenario_id, sim$bin_size)
    out <- lapply(split(sim, key), function(d) {
        tibble::tibble(
            scenario_id = d$scenario_id[1L], bin_size = d$bin_size[1L],
            G = d$G[1L], n = d$n[1L], L = d$L[1L],
            n_samples = d$n_samples[1L], n_bins = d$n_bins[1L],
            n_fits = nrow(d),
            time_median = stats::median(d$wall_time, na.rm = TRUE),
            time_p95 = unname(stats::quantile(d$wall_time, 0.95,
                na.rm = TRUE)),
            time_max = max(d$wall_time, na.rm = TRUE)
        )
    })
    do.call(rbind, unname(out))
}

## =========================================================================
## Tier 1 summaries (A2-2). Rejection is p_value <= alpha; a failed arm (NA
## p-value) counts as a non-rejection in the headline rate (Q6), and the rate
## conditional on success is reported alongside, with the failure rate.
## =========================================================================

ALPHA <- 0.05
ARM_KEY <- c("arm", "test_variant")

.wilson <- function(x, n, z = stats::qnorm(0.975)) {
    p <- x / n
    ctr <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(ctr - hw, ctr + hw)
}

## Newcombe (1998) hybrid-score interval for p1 - p2; `phi` is the paired
## correlation correction (0 for independent samples).
.newcombe_diff <- function(x1, n1, x2, n2, phi = 0) {
    w1 <- .wilson(x1, n1)
    w2 <- .wilson(x2, n2)
    p1 <- x1 / n1
    p2 <- x2 / n2
    d <- p1 - p2
    lo <- d - sqrt((p1 - w1[1])^2 - 2 * phi * (p1 - w1[1]) * (w2[2] - p2) +
        (w2[2] - p2)^2)
    hi <- d + sqrt((w1[2] - p1)^2 - 2 * phi * (w1[2] - p1) * (p2 - w2[1]) +
        (p2 - w2[1])^2)
    c(diff = d, lower = max(lo, -1), upper = min(hi, 1))
}

## Paired rejection indicators (same datasets): McNemar-style comparison.
## Returns discordant counts, the paired difference with a Newcombe paired
## interval, and the exact McNemar p-value.
paired_diff <- function(rej_a, rej_b) {
    N <- length(rej_a)
    a <- sum(rej_a & rej_b)
    b <- sum(rej_a & !rej_b)
    c <- sum(!rej_a & rej_b)
    d <- sum(!rej_a & !rej_b)
    A <- (a + b) * (c + d) * (a + c) * (b + d)
    phi <- 0
    if (A > 0) {
        B <- a * d - b * c
        C <- if (B > N / 2) B - N / 2 else if (B >= 0) 0 else B
        phi <- C / sqrt(A)
    }
    nd <- .newcombe_diff(a + b, N, a + c, N, phi)
    tibble::tibble(n_pairs = N, a_only = b, b_only = c, both = a,
        neither = d, rate_a = (a + b) / N, rate_b = (a + c) / N,
        diff = nd[["diff"]], ci_lower = nd[["lower"]],
        ci_upper = nd[["upper"]],
        mcnemar_p = if (b + c > 0) stats::binom.test(b, b + c, 0.5)$p.value
            else 1)
}

.rej <- function(p) !is.na(p) & p <= ALPHA

## Rate table: one row per group of `by` columns.
rate_table <- function(sim, by) {
    sim$.rej <- .rej(sim$p_value)
    sim$.ok <- sim$status == "ok"
    key <- do.call(paste, c(sim[by], sep = "\r"))
    out <- lapply(split(seq_len(nrow(sim)), key), function(ix) {
        d <- sim[ix, ]
        n <- nrow(d)
        x <- sum(d$.rej)
        ci <- cp_interval(x, n)
        n_ok <- sum(d$.ok)
        ci_c <- cp_interval(sum(d$.rej & d$.ok), n_ok)
        sing <- d$singular[d$.ok]
        r <- tibble::as_tibble(d[1L, by])
        r$n_datasets <- n
        r$n_reject <- x
        r$rate <- x / n
        r$ci_lower <- ci[["lower"]]
        r$ci_upper <- ci[["upper"]]
        r$rate_cond_ok <- if (n_ok) sum(d$.rej & d$.ok) / n_ok else NA_real_
        r$ci_cond_lower <- ci_c[["lower"]]
        r$ci_cond_upper <- ci_c[["upper"]]
        r$failure_rate <- 1 - n_ok / n
        r$singular_rate <- if (length(sing) && !all(is.na(sing)))
            mean(sing, na.rm = TRUE) else NA_real_
        r$kw_perm_random_share <- if (any(!is.na(d$kw_perm_method)))
            mean(d$kw_perm_method == "random_perm", na.rm = TRUE)
            else NA_real_
        r
    })
    do.call(rbind, unname(out))
}

.q7_class <- function(lower, upper) {
    ifelse(upper <= 0.06, "controlled", ifelse(lower > ALPHA, "liberal",
        "neither"))
}

.order_by <- function(tab, cols) tab[do.call(order, tab[cols]), ]

.factorial <- function(sim) sim[sim$family == "factorial", ]

tier1_typeI_table <- function(sim) {
    d <- .factorial(sim)
    d <- d[d$Delta_std == 0, ]
    t <- rate_table(d, c(ARM_KEY, "n", "G"))
    t$q7_class <- .q7_class(t$ci_lower, t$ci_upper)
    .order_by(t, c("arm", "test_variant", "G", "n"))
}

tier1_power_table <- function(sim) {
    d <- .factorial(sim)
    d <- d[d$Delta_std > 0, ]
    .order_by(rate_table(d, c(ARM_KEY, "n", "G", "Delta_std")),
        c("arm", "test_variant", "G", "n", "Delta_std"))
}

## Join two arms' rejection indicators on identical datasets.
.paired_cells <- function(sim, a, b_arms, by_cell) {
    sim$.rej <- .rej(sim$p_value)
    ka <- sim[sim$arm == a$arm & sim$test_variant == a$test_variant, ]
    out <- list()
    for (i in seq_len(nrow(b_arms))) {
        kb <- sim[sim$arm == b_arms$arm[i] &
            sim$test_variant == b_arms$test_variant[i], ]
        m <- merge(ka[c("scenario_id", "bin_size", "replicate_id", ".rej",
            by_cell)], kb[c("scenario_id", "bin_size", "replicate_id",
            ".rej")], by = c("scenario_id", "bin_size", "replicate_id"),
            suffixes = c("_a", "_b"))
        cell <- do.call(paste, c(m[by_cell], sep = "\r"))
        res <- lapply(split(seq_len(nrow(m)), cell), function(ix) {
            r <- tibble::as_tibble(m[ix[1L], by_cell])
            r$comparator_arm <- b_arms$arm[i]
            r$comparator_variant <- b_arms$test_variant[i]
            cbind(r, paired_diff(m$.rej_a[ix], m$.rej_b[ix]))
        })
        out <- c(out, unname(res))
    }
    do.call(rbind, out)
}

## Primary LMM (omnibus) minus each KW variant, positive = LMM rejects more.
tier1_paired_lmm_vs_kw <- function(sim) {
    d <- .factorial(sim)
    b <- data.frame(
        arm = rep(c("kw_chisq", "kw_perm"), each = 3L),
        test_variant = rep(c("omnibus", "dunn_ungated", "dunn_gated"), 2L))
    t <- .paired_cells(d, list(arm = "lmm", test_variant = "omnibus"), b,
        c("n", "G", "Delta_std"))
    t$scenario_kind <- ifelse(t$Delta_std == 0, "typeI", "power")
    names(t)[names(t) == "a_only"] <- "lmm_only_rejects"
    names(t)[names(t) == "b_only"] <- "comparator_only_rejects"
    names(t)[names(t) == "rate_a"] <- "rate_lmm"
    names(t)[names(t) == "rate_b"] <- "rate_comparator"
    .order_by(t, c("comparator_arm", "comparator_variant", "G", "n",
        "Delta_std"))
}

tier1_paired_lmm_vs_twostage <- function(sim) {
    d <- .factorial(sim)
    b <- data.frame(arm = "two_stage", test_variant = "omnibus")
    t <- .paired_cells(d, list(arm = "lmm", test_variant = "omnibus"), b,
        c("n", "G", "Delta_std"))
    t$scenario_kind <- ifelse(t$Delta_std == 0, "typeI", "power")
    names(t)[names(t) == "a_only"] <- "lmm_only_rejects"
    names(t)[names(t) == "b_only"] <- "comparator_only_rejects"
    names(t)[names(t) == "rate_a"] <- "rate_lmm"
    names(t)[names(t) == "rate_b"] <- "rate_comparator"
    .order_by(t, c("G", "n", "Delta_std"))
}

## Each one-at-a-time sweep cell against its reference cell (same n, G = 2,
## Delta_std, reference nuisance). Datasets differ across scenarios (seeds
## derive from the scenario ID), so the difference is an independent-sample
## Newcombe interval, except bin size, which reuses the reference datasets
## and is paired.
tier1_sensitivity_summary <- function(sim) {
    sim$.rej <- .rej(sim$p_value)
    ref <- sim[sim$family == "factorial" & sim$G == 2L, ]
    sw <- sim[sim$family == "sweep", ]
    cells <- unique(sw[c("scenario_id", "bin_size", "sweep_factor",
        "sweep_level", "n", "G", "Delta_std")])
    out <- list()
    for (i in seq_len(nrow(cells))) {
        cl <- cells[i, ]
        rc <- ref[ref$n == cl$n & ref$Delta_std == cl$Delta_std, ]
        vc <- sw[sw$scenario_id == cl$scenario_id &
            sw$bin_size == cl$bin_size, ]
        paired <- cl$sweep_factor == "bin_size"
        for (k in unique(paste(vc$arm, vc$test_variant, sep = "\r"))) {
            kk <- strsplit(k, "\r", fixed = TRUE)[[1L]]
            rv <- vc[vc$arm == kk[1L] & vc$test_variant == kk[2L], ]
            rr <- rc[rc$arm == kk[1L] & rc$test_variant == kk[2L], ]
            if (paired) {
                m <- merge(rr[c("replicate_id", ".rej")],
                    rv[c("replicate_id", ".rej")], by = "replicate_id")
                pd <- paired_diff(m$.rej.y, m$.rej.x)
                nd <- c(pd$diff, pd$ci_lower, pd$ci_upper)
            } else {
                x <- .newcombe_diff(sum(rv$.rej), nrow(rv), sum(rr$.rej),
                    nrow(rr))
                nd <- unname(x)
            }
            out[[length(out) + 1L]] <- tibble::tibble(
                sweep_factor = cl$sweep_factor,
                sweep_level = cl$sweep_level, n = cl$n, G = cl$G,
                Delta_std = cl$Delta_std,
                metric = if (cl$Delta_std == 0) "typeI" else "power",
                arm = kk[1L], test_variant = kk[2L],
                rate_reference = mean(rr$.rej), rate_sweep = mean(rv$.rej),
                diff_vs_reference = nd[1L], ci_lower = nd[2L],
                ci_upper = nd[3L], paired = paired,
                failure_rate_sweep = mean(rv$status != "ok"))
        }
    }
    t <- do.call(rbind, out)
    t$sweep_level_num <- suppressWarnings(as.numeric(t$sweep_level))
    t <- t[order(t$sweep_factor, t$sweep_level_num, t$n, t$Delta_std,
        t$arm, t$test_variant), ]
    t$sweep_level_num <- NULL
    t
}

## Collapse check (PLAN.md 2.5): power of the primary LMM and the two-stage
## arm along the s_ref axis (cv = 1/3) and the cv axis, at n in {3, 5},
## G = 2. `matched_group` is Delta_std / cv, which is Delta / tau. The
## difference is against the reference cell of the same n and matched group
## (s_ref = 6e-6, cv = 1/3, or the cv-sweep cell at cv = 1/3).
tier1_collapse_check <- function(sim) {
    sim$.rej <- .rej(sim$p_value)
    nuis <- sim$arm %in% c("lmm", "two_stage") &
        sim$test_variant == "omnibus" & sim$G == 2L &
        sim$n %in% c(3L, 5L) & sim$bin_size == 5000L &
        sim$sigma_atom == 0.15 & sim$L == 238080 & sim$anchor_frac == 0.5
    axis <- (sim$family %in% c("sweep", "collapse_extra") &
        sim$sweep_factor %in% c("s_ref", "cv")) |
        (sim$family == "factorial" & sim$Delta_std == 1)
    d <- sim[nuis & axis & sim$Delta_std > 0, ]
    d$matched_group <- round(d$Delta_std / d$cv, 3)
    key <- paste(d$arm, d$n, d$s_ref, d$cv, d$Delta_std, sep = "\r")
    out <- lapply(split(seq_len(nrow(d)), key), function(ix) {
        x <- d[ix, ]
        ci <- cp_interval(sum(x$.rej), nrow(x))
        tibble::tibble(arm = x$arm[1L], n = x$n[1L], s_ref = x$s_ref[1L],
            cv = x$cv[1L], Delta_std = x$Delta_std[1L],
            matched_group = x$matched_group[1L], n_datasets = nrow(x),
            power = mean(x$.rej), ci_lower = ci[["lower"]],
            ci_upper = ci[["upper"]])
    })
    t <- do.call(rbind, unname(out))
    refrow <- function(r) {
        t[t$arm == r$arm & t$n == r$n & t$matched_group == r$matched_group &
            t$s_ref == 6e-6 & abs(t$cv - 1 / 3) < 1e-9, ]
    }
    t$diff_vs_reference <- NA_real_
    t$diff_lower <- NA_real_
    t$diff_upper <- NA_real_
    for (i in seq_len(nrow(t))) {
        rf <- refrow(t[i, ])
        if (nrow(rf) == 1L && !(t$s_ref[i] == 6e-6 &&
            abs(t$cv[i] - 1 / 3) < 1e-9)) {
            nd <- .newcombe_diff(round(t$power[i] * t$n_datasets[i]),
                t$n_datasets[i], round(rf$power * rf$n_datasets),
                rf$n_datasets)
            t$diff_vs_reference[i] <- nd[["diff"]]
            t$diff_lower[i] <- nd[["lower"]]
            t$diff_upper[i] <- nd[["upper"]]
        }
    }
    t$agrees_with_reference <- ifelse(is.na(t$diff_lower), NA,
        t$diff_lower <= 0 & t$diff_upper >= 0)
    t[order(t$arm, t$n, t$matched_group, t$s_ref, t$cv), ]
}

## M3: Type I of the primary LMM by anchor position (Delta_std = 0, G = 2).
tier1_m3_anchor <- function(sim) {
    d <- sim[sim$arm == "lmm" & sim$Delta_std == 0 & sim$G == 2L &
        sim$n %in% c(3L, 5L) & sim$bin_size == 5000L &
        sim$cv == 1 / 3 & sim$sigma_atom == 0.15 & sim$L == 238080 &
        sim$s_ref == 6e-6 &
        (sim$family == "factorial" | sim$sweep_factor == "anchor_frac"), ]
    t <- rate_table(d, c("anchor_frac", "n"))
    pool <- rate_table(transform(d, n = 0L), c("anchor_frac", "n"))
    t <- rbind(t, pool)
    t$q7_class <- .q7_class(t$ci_lower, t$ci_upper)
    ## Heterogeneity across anchors (chi-square test of homogeneity, exact
    ## by simulation) within each n; n = 0 pools n = 3 and 5.
    t$homogeneity_p <- vapply(t$n, function(nn) {
        s <- t[t$n == nn, ]
        tryCatch(stats::fisher.test(cbind(s$n_reject,
            s$n_datasets - s$n_reject), simulate.p.value = TRUE,
            B = 1e5)$p.value, error = function(e) NA_real_)
    }, numeric(1L))
    t[order(t$n, t$anchor_frac), ]
}

## Failure counts and timing per arm and design size.
tier1_failures_timing <- function(sim) {
    sim$.ok <- sim$status == "ok"
    key <- paste(sim$arm, sim$test_variant, sim$bin_size, sep = "\r")
    out <- lapply(split(seq_len(nrow(sim)), key), function(ix) {
        d <- sim[ix, ]
        tibble::tibble(arm = d$arm[1L], test_variant = d$test_variant[1L],
            bin_size = d$bin_size[1L], n_rows = nrow(d),
            n_failed = sum(!d$.ok), failure_rate = mean(!d$.ok),
            first_error = if (any(!d$.ok)) d$status[!d$.ok][1L]
                else NA_character_,
            wall_median_s = stats::median(d$wall_time, na.rm = TRUE),
            wall_p95_s = unname(stats::quantile(d$wall_time, 0.95,
                na.rm = TRUE)),
            wall_total_s = sum(d$wall_time, na.rm = TRUE))
    })
    do.call(rbind, unname(out))
}

write_tier1_tables <- function(sim, out_dir = "results/explore/tier1") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    tabs <- list(
        "tier1-typeI-by-arm.csv" = tier1_typeI_table(sim),
        "tier1-power-by-arm.csv" = tier1_power_table(sim),
        "tier1-paired-lmm-vs-kw.csv" = tier1_paired_lmm_vs_kw(sim),
        "tier1-paired-lmm-vs-twostage.csv" = tier1_paired_lmm_vs_twostage(sim),
        "tier1-sensitivity-summary.csv" = tier1_sensitivity_summary(sim),
        "tier1-collapse-check.csv" = tier1_collapse_check(sim),
        "tier1-m3-anchor.csv" = tier1_m3_anchor(sim),
        "tier1-failures-timing.csv" = tier1_failures_timing(sim)
    )
    paths <- file.path(out_dir, names(tabs))
    for (i in seq_along(tabs)) {
        utils::write.csv(tabs[[i]], paths[i], row.names = FALSE)
    }
    rds <- file.path(out_dir, "tier1-replicates.rds")
    saveRDS(sim, rds, compress = "xz")
    c(paths, rds)
}
