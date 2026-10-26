## =========================================================================
## Pipeline for the simulation study (PLAN.md 5.1, 5.5). Run from
## simulations/.
##
##   CDA_SIM_TIER     tier to build: "pilot", "explore" (Tier 1), "a23" (A2-3; use a
##                   separate store: tar_make_future(store = "_targets_a23")) or
##                   the A2-4 stages "a24typeI", "a24power", "a24sens" (one store
##                   each, `_targets_a24typeI` etc.; tables are written by
##                   scripts/a24_driver.R, which also runs the stage checks)
##   CDA_SIM_R        replicates per scenario (default 20)
##   CDA_SIM_BATCH    replicates per batch (default 20)
##   CDA_SIM_WORKERS  1 (default): targets::tar_make(), serial.
##                    >1: future multisession; run targets::tar_make_future().
##
## The unit of work is a batch (one scenario x bin size x block of
## replicates); dynamic branching goes over batches. Only the Tier 0 targets
## are declared here: scenario grid, generator and all five analysis arms (run
## together inside `sim_batch`, so every arm is paired on the same dataset),
## the timing summary, and the Tier 1 summary tables (written when the tier is
## "explore").
## =========================================================================

library(targets)

source("R/setup.R")
load_sim()

tier    <- Sys.getenv("CDA_SIM_TIER", "pilot")
R       <- as.integer(Sys.getenv("CDA_SIM_R", "20"))
batch   <- as.integer(Sys.getenv("CDA_SIM_BATCH", "20"))
workers <- as.integer(Sys.getenv("CDA_SIM_WORKERS", "1"))

if (workers > 1L) {
    ## Parallelism is at the branch level only; pin BLAS/OpenMP threads so
    ## workers do not oversubscribe the machine.
    ## RhpcBLASctl is not in the renv library, so the threads are pinned by
    ## environment variable, before any worker starts.
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
    future::plan(future::multisession, workers = workers)
}

tar_option_set(
    packages = c("tibble", "data.table"),
    error = "continue"
)

list(
    ## Files that define the statistics. A change to any invalidates the
    ## downstream simulation branches; unrelated package edits do not.
    tar_target(
        pkg_methods,
        file.path("..", "R", c("decay_tests.R", "null_distribution.R",
            "decay_slopes.R", "coverage_io.R")),
        format = "file"
    ),
    tar_target(grid, build_grid(tier)),
    tar_target(batches, expand_batches(grid, R, batch)),
    tar_target(
        sim_batch,
        if (grepl("^a24", tier)) {
            run_batch_slim(batches, pkg_methods, arms = a24_arms_for(batches))
        }
        else run_batch(batches, pkg_methods),
        pattern = map(batches)
    ),
    tar_target(sim_all, sim_batch),
    if (!grepl("^a24", tier)) {
        tar_target(timing_summary, summarize_timing(sim_all))
    },
    if (tier == "a23") {
        tar_target(
            a23_tables,
            write_a23_tables(sim_all, file.path("results", "explore", "a2-3")),
            format = "file"
        )
    },
    if (tier == "explore") {
        tar_target(
            tier1_tables,
            write_tier1_tables(sim_all,
                file.path("results", "explore", "tier1")),
            format = "file"
        )
    }
)
