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
    ## Numeric, not integer: the product of the four margins overflows an
    ## integer at N around 2,000 and up.
    a <- as.numeric(sum(rej_a & rej_b))
    b <- as.numeric(sum(rej_a & !rej_b))
    c <- as.numeric(sum(!rej_a & rej_b))
    d <- as.numeric(sum(!rej_a & !rej_b))
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

## =========================================================================
## A2-3 tables.
## =========================================================================

.bind_fill <- function(dfs) {
    cols <- unique(unlist(lapply(dfs, names)))
    dfs <- lapply(dfs, function(d) {
        for (cn in setdiff(cols, names(d))) d[[cn]] <- NA
        d[cols]
    })
    do.call(rbind, dfs)
}

.in_family <- function(sim, fam) {
    vapply(strsplit(sim$family, ";", fixed = TRUE), function(x) fam %in% x,
        logical(1L))
}

## Anchor 0.25 cell. The 200 replicates shared with Tier 1 are where the cell
## was flagged (selected because it looked high), so the fresh replicates
## 201 to R are reported separately as the unselected confirmation.
a23_anchor_table <- function(sim) {
    d <- sim[.in_family(sim, "anchor_025") & sim$bin_size == 5000L, ]
    sub <- list(all = d, fresh_201_on = d[d$replicate_id > 200L, ],
        tier1_1_200 = d[d$replicate_id <= 200L, ])
    out <- lapply(names(sub), function(nm) {
        x <- sub[[nm]]
        x$subset <- nm
        t <- rate_table(x, c("subset", ARM_KEY))
        t$q7_class <- .q7_class(t$ci_lower, t$ci_upper)
        t
    })
    do.call(rbind, out)
}

## 1 kb vs the Tier 1 5 kb reference cells (joined by scenario_id and
## replicate_id). Tier 1 has replicates 1 to 200 only, so the paired
## comparison uses those; the unpaired comparison uses all 1 kb replicates.
## `extra_5kb` (optional): further 5 kb replicates of the same reference
## scenarios (the A2-3 supplement, replicates 201 to 1000), appended to the
## Tier 1 rows. Only arm/variants present in Tier 1 are compared then, so
## every row of the table rests on the same pairs.
a23_1kb_paired <- function(sim, tier1, extra_5kb = NULL) {
    k1 <- sim[.in_family(sim, "one_kb") & sim$bin_size == 1000L, ]
    t5 <- tier1[tier1$bin_size == 5000L & tier1$family == "factorial", ]
    if (!is.null(extra_5kb)) {
        cols <- c("scenario_id", "replicate_id", "bin_size", "arm",
            "test_variant", "p_value", "status", "Delta_std")
        t5 <- rbind(t5[cols], extra_5kb[extra_5kb$bin_size == 5000L, cols])
        stopifnot(!anyDuplicated(t5[c("scenario_id", "replicate_id", "arm",
            "test_variant")]))
        k1 <- k1[paste(k1$arm, k1$test_variant) %in%
            paste(tier1$arm, tier1$test_variant), ]
    }
    k1$.rej <- .rej(k1$p_value)
    t5$.rej <- .rej(t5$p_value)
    out <- list()
    cells <- unique(k1[c("scenario_id", "Delta_std")])
    for (i in seq_len(nrow(cells))) {
        sid <- cells$scenario_id[i]
        for (k in unique(paste(k1$arm[k1$scenario_id == sid],
            k1$test_variant[k1$scenario_id == sid], sep = "\r"))) {
            kk <- strsplit(k, "\r", fixed = TRUE)[[1L]]
            a <- k1[k1$scenario_id == sid & k1$arm == kk[1L] &
                k1$test_variant == kk[2L], ]
            b <- t5[t5$scenario_id == sid & t5$arm == kk[1L] &
                t5$test_variant == kk[2L], ]
            if (!nrow(b)) next
            m <- merge(a[c("replicate_id", ".rej")],
                b[c("replicate_id", ".rej")], by = "replicate_id")
            pd <- paired_diff(m$.rej.x, m$.rej.y)
            un <- .newcombe_diff(sum(a$.rej), nrow(a), sum(b$.rej), nrow(b))
            out[[length(out) + 1L]] <- tibble::tibble(
                Delta_std = cells$Delta_std[i],
                metric = if (cells$Delta_std[i] == 0) "typeI" else "power",
                arm = kk[1L], test_variant = kk[2L],
                n_1kb = nrow(a), rate_1kb = mean(a$.rej),
                n_5kb_tier1 = nrow(b), rate_5kb_tier1 = mean(b$.rej),
                n_pairs = pd$n_pairs,
                rate_1kb_on_pairs = pd$rate_a, rate_5kb_on_pairs = pd$rate_b,
                paired_diff_1kb_minus_5kb = pd$diff,
                paired_ci_lower = pd$ci_lower, paired_ci_upper = pd$ci_upper,
                discordant_1kb_only = pd$a_only,
                discordant_5kb_only = pd$b_only, mcnemar_p = pd$mcnemar_p,
                unpaired_diff_all = un[["diff"]],
                unpaired_ci_lower = un[["lower"]],
                unpaired_ci_upper = un[["upper"]],
                failure_rate_1kb = mean(a$status != "ok"))
        }
    }
    do.call(rbind, out)
}

## All arms at the single slope-only reference cell (G = 2, n = 5, Delta_std 0).
a23_slope_only_table <- function(sim) {
    d <- sim[.in_family(sim, "slope_only") & sim$bin_size == 5000L, ]
    t <- rate_table(d, ARM_KEY)
    t$q7_class <- .q7_class(t$ci_lower, t$ci_upper)
    t
}

## Pooled-variance vs Welch two-stage, and pooled vs the primary LMM, paired
## on the same datasets in every validation cell.
a23_pooled_vs_welch <- function(sim) {
    d <- sim[(.in_family(sim, "pooled") | .in_family(sim, "pooled_extra")) &
        sim$bin_size == 5000L & sim$test_variant == "omnibus", ]
    d$.rej <- .rej(d$p_value)
    out <- list()
    for (sid in unique(d$scenario_id)) {
        x <- d[d$scenario_id == sid, ]
        p <- x[x$arm == "two_stage_pooled", ]
        for (cmp in c("two_stage", "lmm")) {
            o <- x[x$arm == cmp, ]
            m <- merge(p[c("replicate_id", ".rej", "p_value")],
                o[c("replicate_id", ".rej", "p_value")], by = "replicate_id")
            pd <- paired_diff(m$.rej.x, m$.rej.y)
            out[[length(out) + 1L]] <- cbind(
                tibble::tibble(G = p$G[1L], n = p$n[1L],
                    Delta_std = p$Delta_std[1L],
                    metric = if (p$Delta_std[1L] == 0) "typeI" else "power",
                    comparator = cmp,
                    n_identical_p = sum(abs(m$p_value.x - m$p_value.y) < 1e-8),
                    max_abs_p_diff = max(abs(m$p_value.x - m$p_value.y))),
                pd)
        }
    }
    names(out) <- NULL
    t <- do.call(rbind, out)
    names(t)[names(t) == "a_only"] <- "pooled_only_rejects"
    names(t)[names(t) == "b_only"] <- "comparator_only_rejects"
    names(t)[names(t) == "rate_a"] <- "rate_pooled"
    names(t)[names(t) == "rate_b"] <- "rate_comparator"
    t[order(t$comparator, t$G, t$n, t$Delta_std), ]
}

## Permutation-Dunn validation at G = 3, n = 5 (Delta_std 0 and 1): finite
## p-values, method counts, rejection rates, and agreement with the
## normal-approximation Dunn rows on the same datasets.
a23_dunn_perm_validation <- function(sim) {
    d <- sim[.in_family(sim, "dunn_perm") & sim$bin_size == 5000L &
        sim$arm %in% c("kw_chisq", "kw_perm"), ]
    d$.rej <- .rej(d$p_value)
    rates <- list()
    for (sid in unique(d$scenario_id)) {
        for (k in unique(paste(d$arm[d$scenario_id == sid],
            d$test_variant[d$scenario_id == sid], sep = "\r"))) {
            kk <- strsplit(k, "\r", fixed = TRUE)[[1L]]
            x <- d[d$scenario_id == sid & d$arm == kk[1L] &
                d$test_variant == kk[2L], ]
            ci <- cp_interval(sum(x$.rej), nrow(x))
            rates[[length(rates) + 1L]] <- tibble::tibble(
                row_type = "rate", G = x$G[1L], n = x$n[1L],
                Delta_std = x$Delta_std[1L], arm = kk[1L],
                test_variant = kk[2L], n_rows = nrow(x),
                n_nonfinite_p = sum(!is.finite(x$p_value)),
                n_reject = sum(x$.rej), rate = mean(x$.rej),
                ci_lower = ci[["lower"]], ci_upper = ci[["upper"]],
                n_method_exact_enum = sum(x$dunn_perm_method == "exact_enum",
                    na.rm = TRUE),
                n_method_random_perm = sum(x$dunn_perm_method == "random_perm",
                    na.rm = TRUE))
        }
    }
    pairs <- list()
    pk <- d[d$arm == "kw_perm", ]
    for (sid in unique(pk$scenario_id)) {
        x <- pk[pk$scenario_id == sid, ]
        for (cmp in list(c("dunn_perm_ungated", "dunn_ungated"),
            c("dunn_perm_gated", "dunn_gated"))) {
            a <- x[x$test_variant == cmp[1L], ]
            b <- x[x$test_variant == cmp[2L], ]
            m <- merge(a[c("replicate_id", ".rej")],
                b[c("replicate_id", ".rej")], by = "replicate_id")
            pd <- paired_diff(m$.rej.x, m$.rej.y)
            pairs[[length(pairs) + 1L]] <- cbind(
                tibble::tibble(row_type = "paired_perm_vs_approx",
                    G = x$G[1L], n = x$n[1L], Delta_std = x$Delta_std[1L],
                    arm = "kw_perm", test_variant = cmp[1L],
                    comparator_variant = cmp[2L]), pd)
        }
    }
    .bind_fill(c(rates, pairs))
}

## Writes the A2-3 tables. The 1 kb paired table joins the supplemental 5 kb
## replicates 201 to 1000 (`extra_5kb_path`, from scripts/a23_5kb_supplement.R)
## to Tier 1 for the full R = 1000 pairing; a missing file is an error, not a
## silent fall back to the R = 200 pairing.
write_a23_tables <- function(sim, out_dir = "results/explore/a2-3",
                             tier1_path = "results/explore/tier1/tier1-replicates.rds",
                             extra_5kb_path = file.path(out_dir,
                                 "a2-3-5kb-supplement-replicates.rds")) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    tier1 <- readRDS(tier1_path)
    if (!file.exists(extra_5kb_path)) {
        stop("Supplemental 5 kb replicates not found: ", extra_5kb_path,
            ". Run scripts/a23_5kb_supplement.R first.")
    }
    extra_5kb <- readRDS(extra_5kb_path)
    tabs <- list(
        "a2-3-anchor-025-typeI.csv" = a23_anchor_table(sim),
        "a2-3-1kb-paired.csv" = a23_1kb_paired(sim, tier1, extra_5kb),
        "a2-3-slope-only-high-R.csv" = a23_slope_only_table(sim),
        "a2-3-pooled-vs-welch.csv" = a23_pooled_vs_welch(sim),
        "a2-3-dunn-perm-validation.csv" = a23_dunn_perm_validation(sim),
        "a2-3-failures-timing.csv" = tier1_failures_timing(sim)
    )
    paths <- file.path(out_dir, names(tabs))
    for (i in seq_along(tabs)) {
        utils::write.csv(tabs[[i]], paths[i], row.names = FALSE)
    }
    caption <- file.path(out_dir, "a2-3-1kb-paired-CAPTION.txt")
    writeLines(a23_1kb_paired_caption(), caption)
    rds <- file.path(out_dir, "a2-3-replicates.rds")
    saveRDS(sim, rds, compress = "xz")
    c(paths, caption, rds)
}

## Caption for a2-3-1kb-paired.csv, written beside it (a CSV cannot carry
## comment lines without breaking read.csv).
a23_1kb_paired_caption <- function() {
    c("a2-3-1kb-paired.csv: 1 kb vs 5 kb bins, paired on the same datasets.",
        "",
        "Cells: G = 2, n = 5, reference nuisance, Delta_std 0 (Type I), 0.5, 1, 1.5 (power).",
        "Pairing: 1 kb replicates 1 to 1000 joined by scenario_id and replicate_id to",
        "5 kb replicates 1 to 200 (Tier 1) plus 201 to 1000 (A2-3 supplement,",
        "scripts/a23_5kb_supplement.R). n_pairs = 1000 in every row.",
        "Columns: paired_diff = 1 kb minus 5 kb rejection rate; paired_ci_* = Newcombe",
        "paired-difference 95% interval; mcnemar_p = exact McNemar test; unpaired_*",
        "compare the two marginal rates.",
        "",
        "Arm variants included (the nine present in Tier 1): lmm; lmm_slope_only_CAUTIONARY;",
        "kw_chisq (omnibus, dunn_ungated, dunn_gated); kw_perm (omnibus, dunn_ungated,",
        "dunn_gated); two_stage (omnibus).",
        "Deferred to Tier 2: two_stage_pooled and kw_perm dunn_perm_ungated/dunn_perm_gated.",
        "They exist only at 5 kb replicates 201 to 1000 (Tier 1 predates them), so pairing",
        "them would rest on 800, not 1000, pairs and would not match the rows above.")
}
