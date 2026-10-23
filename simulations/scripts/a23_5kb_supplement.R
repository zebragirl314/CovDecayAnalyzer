## A2-3 supplement. Run from simulations/:  Rscript scripts/a23_5kb_supplement.R
## 5 kb reference cells (G = 2, n = 5, Delta_std {0, .5, 1, 1.5}) at replicates
## 201 to 1000, to pair with the A2-3 1 kb R = 1000 run. Tier 1 holds 5 kb
## replicates 1 to 200. Seeds derive from (scenario_id, replicate_id), so the
## new replicates are the same datasets the 1 kb run analyzed.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("R/setup.R")
load_sim()
out_dir <- "results/explore/a2-3"
workers <- 10L

grid <- do.call(rbind, lapply(c(0, 0.5, 1, 1.5), function(d) {
    scenario_to_row(make_scenario(G = 2L, n = 5L, Delta_std = d,
        bin_size = 5000L))
}))
grid$family <- "factorial"
tier1 <- readRDS("results/explore/tier1/tier1-replicates.rds")
a23 <- readRDS(file.path(out_dir, "a2-3-replicates.rds"))
stopifnot(all(grid$scenario_id %in% tier1$scenario_id),
    all(grid$scenario_id %in% a23$scenario_id[a23$bin_size == 1000L]))

## Reproducibility check: current code on replicates 1 to 10 against Tier 1.
chk <- grid
chk$rep_start <- 1L
chk$rep_end <- 10L
chk_res <- do.call(rbind, parallel::mclapply(seq_len(nrow(chk)),
    function(i) run_batch(chk[i, ]), mc.cores = 4L))
m <- merge(chk_res[c("scenario_id", "replicate_id", "arm", "test_variant",
    "p_value")], tier1[tier1$bin_size == 5000L,
    c("scenario_id", "replicate_id", "arm", "test_variant", "p_value")],
    by = c("scenario_id", "replicate_id", "arm", "test_variant"))
cat(sprintf("Repro check: %d rows matched to Tier 1, max |dp| = %.3g\n",
    nrow(m), max(abs(m$p_value.x - m$p_value.y), na.rm = TRUE)))
print(table(m$arm, m$test_variant))

b <- expand_batches(grid, 1000L, 50L)
b <- b[b$rep_start >= 201L, ]
cat(nrow(b), "batches\n")
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(b)), function(i) run_batch(b[i, ]),
    mc.cores = workers, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1L), "try-error")
stopifnot(!any(bad))
sim5 <- do.call(rbind, res)
cat("Elapsed (min):", format(difftime(Sys.time(), t0, units = "mins")), "\n")
saveRDS(sim5, file.path(out_dir, "a2-3-5kb-supplement-replicates.rds"),
    compress = "xz")

paired <- a23_1kb_paired(a23, tier1, extra_5kb = sim5)
utils::write.csv(paired, file.path(out_dir, "a2-3-1kb-paired.csv"),
    row.names = FALSE)
cat("failed rows:", sum(sim5$status != "ok"), "of", nrow(sim5), "\n")
