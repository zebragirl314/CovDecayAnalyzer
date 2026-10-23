## Supplement to the Tier 0 pilot: wall time of the primary test at 1 kb
## bins across genome lengths and designs. The main pilot timed only one 1 kb
## design (the largest, about 57 s per fit), which is not enough to cost the
## paired 1 kb sensitivity (PLAN.md 4.1). R = 3 fits per cell, serial.
## Run from simulations/:  Rscript scripts/pilot_timing_1kb.R

source("R/setup.R")
load_sim()

out_dir <- Sys.getenv("CDA_PILOT_OUT", "results/explore/pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
R <- as.integer(Sys.getenv("CDA_PILOT_R_1KB", "3"))

cells <- expand.grid(L = c(50000, 100000, 238080, 500000),
    design = c("G2n3", "G2n5", "G4n10"), stringsAsFactors = FALSE)
cells$G <- ifelse(cells$design == "G4n10", 4L, 2L)
cells$n <- c(G2n3 = 3L, G2n5 = 5L, G4n10 = 10L)[cells$design]
## The 500 kb, G = 4, n = 10 design is already timed in the main pilot.
cells <- cells[!(cells$L == 500000 & cells$design == "G4n10"), ]

invisible(analyze_lmm(generate_replicate(make_scenario(), 1L)))  # warm-up
rows <- lapply(seq_len(nrow(cells)), function(i) {
    sc <- make_scenario(L = cells$L[i], G = cells$G[i], n = cells$n[i],
        bin_size = 1000L)
    res <- do.call(rbind, lapply(seq_len(R), function(r) {
        analyze_lmm(generate_replicate(sc, r))
    }))
    cat(sprintf("L=%6.0f G=%d n=%2d: %d obs, median %.2f s\n", sc$L, sc$G,
        cells$n[i], res$n_bins[1L], stats::median(res$wall_time)))
    tibble::tibble(L = sc$L, G = sc$G, n = cells$n[i],
        n_samples = res$n_samples[1L], bins_per_sample = sc$L / 1000,
        obs = res$n_bins[1L], fits = R,
        median_s = stats::median(res$wall_time),
        max_s = max(res$wall_time), n_failed = sum(res$status != "ok"))
})
out <- do.call(rbind, rows)
utils::write.csv(out, file.path(out_dir, "timing-1kb-supplement.csv"),
    row.names = FALSE)
