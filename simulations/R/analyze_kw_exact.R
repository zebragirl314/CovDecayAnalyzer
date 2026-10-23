## KW on per-sample slopes with an exact permutation p-value, plus
## Dunn-Holm ungated and gated. Returns 3 rows: arm = "kw_exact",
## "kw_exact_dunn_ungated", "kw_exact_dunn_gated". The Dunn step uses the
## same normal approximation as the chi-square arm; only the omnibus
## p-value (and so the gate) is exact.
##
## The permutation distribution enumerates every assignment of the ranks to
## groups when there are at most KW_EXACT_MAX assignments (the multinomial
## coefficient). Beyond that it falls back to KW_EXACT_B random
## permutations with p = (k + 1) / (B + 1), which is NOT exact; with
## balanced groups that happens for G = 2, n >= 9 and G = 3, n >= 5.
KW_EXACT_MAX <- 2e5
KW_EXACT_B <- 2e4

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

kw_exact_p <- function(slope, group) {
    N <- length(slope)
    ng <- as.integer(table(group))
    rk <- rank(slope)
    obs_sums <- as.numeric(tapply(rk, group, sum))
    H <- .kw_stat(obs_sums, ng, N)
    n_assign <- exp(lgamma(N + 1) - sum(lgamma(ng + 1)))
    if (n_assign <= KW_EXACT_MAX) {
        sums <- .enum_rank_sums(rk, ng)
        Hp <- apply(sums, 1L, .kw_stat, ng = ng, N = N)
        p <- mean(Hp >= H - 1e-9)
    } else {
        g <- as.integer(group)
        Hp <- replicate(KW_EXACT_B, {
            s <- as.numeric(tapply(rk, sample(g), sum))
            .kw_stat(s, ng, N)
        })
        p <- (sum(Hp >= H - 1e-9) + 1) / (KW_EXACT_B + 1)
    }
    list(p = p, H = H)
}

analyze_kw_exact <- function(data,
                             scenario_id = attr(data, "scenario_id"),
                             replicate_id = attr(data, "replicate_id"),
                             bin_size = attr(data, "bin_size")) {
    ensure_package()
    .kw_rows("kw_exact", data, kw_exact_p, scenario_id, replicate_id,
        bin_size)
}
