## KW on per-sample slopes with a permutation p-value, plus Dunn-Holm
## post-hocs. Returns 5 rows: arm = "kw_perm" with test_variant omnibus,
## dunn_ungated, dunn_gated, dunn_perm_ungated, dunn_perm_gated (internal arm
## strings "kw_perm", "kw_perm_dunn_ungated", "..._dunn_gated",
## "..._dunn_perm_ungated", "..._dunn_perm_gated"). The dunn_ungated/gated
## rows use the normal approximation (as in the chi-square arm; only their
## gate is exact). The dunn_perm rows permute the pairwise Dunn |z| over the
## same permutations as the omnibus (so they use the same enumeration or
## random draws, recorded in `dunn_perm_method`), with Holm adjustment; the
## gated row uses the permutation omnibus p < 0.05. For G = 2 the permutation
## Dunn p equals the omnibus permutation p.
##
## The permutation distribution enumerates every assignment of the ranks to
## groups when there are at most KW_PERM_MAX assignments (the multinomial
## coefficient): method "exact_enum". Beyond that it falls back to KW_PERM_B
## random permutations with p = (k + 1) / (B + 1), which is NOT exact: method
## "random_perm". With balanced groups that happens for G = 2, n >= 11;
## G = 3, n >= 5; and G = 4, n >= 3. The method is stored per fit in the
## `kw_perm_method` column. Random permutations run under a seed derived from
## (scenario_id, replicate_id), so a fit is reproducible and the caller's RNG
## stream is untouched.
KW_PERM_MAX <- 2e5
KW_PERM_B <- 2e4

.kw_stat <- function(rank_sums, ng, N) {
    12 / (N * (N + 1)) * sum(rank_sums^2 / ng) - 3 * (N + 1)
}

## Rank sums of every assignment of `ranks` to groups of sizes `ng`.
.enum_rank_sums <- function(ranks, ng) {
    if (length(ng) == 1L) {
        return(matrix(sum(ranks), 1L, 1L))
    }
    ch <- utils::combn(length(ranks), ng[1L])
    do.call(rbind, lapply(seq_len(ncol(ch)), function(j) {
        rest <- .enum_rank_sums(ranks[-ch[, j]], ng[-1L])
        cbind(sum(ranks[ch[, j]]), rest)
    }))
}

kw_perm_p <- function(slope, group, seed = NULL) {
    N <- length(slope)
    ng <- as.integer(table(group))
    rk <- rank(slope)
    obs_sums <- as.numeric(tapply(rk, group, sum))
    H <- .kw_stat(obs_sums, ng, N)
    n_assign <- exp(lgamma(N + 1) - sum(lgamma(ng + 1)))
    ## `sums`: one row per permutation, one column per group (rank sums).
    if (n_assign <= KW_PERM_MAX) {
        sums <- .enum_rank_sums(rk, ng)
        method <- "exact_enum"
    } else {
        g <- as.integer(group)
        draw <- function() {
            t(replicate(KW_PERM_B, as.numeric(tapply(rk, sample(g), sum))))
        }
        sums <- if (is.null(seed)) draw() else with_seed(seed, draw())
        method <- "random_perm"
    }
    tail_p <- function(stat_perm, stat_obs) {
        k <- sum(stat_perm >= stat_obs - 1e-9)
        if (method == "exact_enum") k / length(stat_perm)
        else (k + 1) / (nrow(sums) + 1)
    }
    Hp <- apply(sums, 1L, .kw_stat, ng = ng, N = N)
    p <- tail_p(Hp, H)
    ## Permutation Dunn: the same permutations, |z| for every pair from the
    ## rank-sum differences (pooled-rank normal-approximation SE, which is
    ## constant across permutations for fixed group sizes), Holm adjusted.
    G <- length(ng)
    pr <- utils::combn(G, 2L)
    se <- sqrt(N * (N + 1) / 12 * (1 / ng[pr[1L, ]] + 1 / ng[pr[2L, ]]))
    rbar <- function(sm) sweep(sm, 2L, ng, "/")
    zmat <- function(sm) {
        rb <- rbar(sm)
        abs(rb[, pr[1L, ], drop = FALSE] - rb[, pr[2L, ], drop = FALSE]) /
            rep(se, each = nrow(sm))
    }
    zp <- zmat(sums)
    zo <- drop(zmat(matrix(obs_sums, 1L)))
    raw <- vapply(seq_len(ncol(pr)), function(j) tail_p(zp[, j], zo[j]),
        numeric(1L))
    list(p = p, H = H, method = method,
        dunn_perm = stats::p.adjust(raw, method = "holm"))
}

analyze_kw_perm <- function(data,
                            scenario_id = attr(data, "scenario_id"),
                            replicate_id = attr(data, "replicate_id"),
                            bin_size = attr(data, "bin_size")) {
    ensure_package()
    seed <- derive_seed(scenario_id %||% "none", replicate_id %||% 0L)
    kw_fun <- function(slope, group) kw_perm_p(slope, group, seed = seed)
    .kw_rows("kw_perm", data, kw_fun, scenario_id, replicate_id, bin_size,
        perm_dunn = TRUE)
}
