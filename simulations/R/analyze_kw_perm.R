## KW on per-sample slopes with a permutation p-value, plus Dunn-Holm
## ungated and gated. Returns 3 rows: arm = "kw_perm", "kw_perm_dunn_ungated",
## "kw_perm_dunn_gated". The Dunn step uses the same normal approximation as
## the chi-square arm; only the omnibus p-value (and so the gate) uses the
## permutation distribution.
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
    if (n_assign <= KW_PERM_MAX) {
        sums <- .enum_rank_sums(rk, ng)
        Hp <- apply(sums, 1L, .kw_stat, ng = ng, N = N)
        p <- mean(Hp >= H - 1e-9)
        method <- "exact_enum"
    } else {
        g <- as.integer(group)
        draw <- function() {
            replicate(KW_PERM_B, {
                s <- as.numeric(tapply(rk, sample(g), sum))
                .kw_stat(s, ng, N)
            })
        }
        Hp <- if (is.null(seed)) draw() else with_seed(seed, draw())
        p <- (sum(Hp >= H - 1e-9) + 1) / (KW_PERM_B + 1)
        method <- "random_perm"
    }
    list(p = p, H = H, method = method)
}

analyze_kw_perm <- function(data,
                            scenario_id = attr(data, "scenario_id"),
                            replicate_id = attr(data, "replicate_id"),
                            bin_size = attr(data, "bin_size")) {
    ensure_package()
    seed <- derive_seed(scenario_id %||% "none", replicate_id %||% 0L)
    kw_fun <- function(slope, group) kw_perm_p(slope, group, seed = seed)
    .kw_rows("kw_perm", data, kw_fun, scenario_id, replicate_id, bin_size)
}
