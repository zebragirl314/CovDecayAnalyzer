## A2-4 driver. Run from simulations/ under caffeinate:
##   caffeinate -dims Rscript scripts/a24_driver.R [typeI power sensitivities]
##
## For each stage: run the targets pipeline (future multisession, 10
## workers, BLAS pinned) as a child process; sample worker memory and cache
## size every 30 s and write the 5-minute resource trace; stop the child on a
## panic exit; then write the stage tables, run `check_stage_results()`, append
## a status block to the running session log, and continue only if the check
## says so. A hard stop or panic writes the interrupted summary stub and exits.
suppressMessages({
    library(data.table)
    library(targets)
})
source("R/setup.R")
load_sim()

DATE <- "2026-10-23"
LOG <- file.path("..", "session-logs", sprintf("session-%s-a2-4.md", DATE))
INTERRUPTED <- file.path("..", "session-logs",
    sprintf("session-%s-a2-4-interrupted.md", DATE))
TRACE <- file.path("..", "session-logs", "a2-4-resource-trace.csv")
WORKERS <- 10L
SESSION_CAP_H <- 24L  # total session wall-time panic threshold (was 20)
BATCH <- 100L

STAGES <- list(
    typeI = list(tier = "a24typeI", budget_h = 6.5, cap_h = 10,
        out = "results/tier2/typeI", writer = write_a24_typeI),
    power = list(tier = "a24power", budget_h = 10, cap_h = 15,
        out = "results/tier2/power", writer = write_a24_power),
    sensitivities = list(tier = "a24sens", budget_h = 6, cap_h = 9,
        out = "results/tier2/sensitivities", writer = write_a24_sens),
    ## S3 follow-up at sigma_atom {0.3, 0.5}; checked as "sensitivities".
    s3rerun = list(tier = "a24s3rerun", budget_h = 1, cap_h = 1.5,
        out = "results/tier2/sensitivities", writer = write_a24_s3rerun,
        check = "sensitivities", log = "session-2026-10-24-s3-rerun"),
    ## Controlled baseline (sigma_ratio = 1) for the S3 rerun.
    s3baseline = list(tier = "a24s3baseline", budget_h = 0.5, cap_h = 0.75,
        out = "results/tier2/sensitivities", writer = write_a24_s3baseline,
        check = "sensitivities", log = "session-2026-10-24-s3-rerun-baseline")
)

args <- commandArgs(trailingOnly = TRUE)
todo <- if (length(args)) args else names(STAGES)
session_start <- Sys.time()

say <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n",
    sep = " ", file = stderr())

append_log <- function(text) {
    cat(text, file = LOG, append = TRUE, sep = "")
}

## Worker RSS (MB, max over workers) and number of workers seen.
sample_workers <- function() {
    ps <- suppressWarnings(system2("ps", c("-axo", "rss=,command="),
        stdout = TRUE, stderr = FALSE))
    w <- ps[grepl("workRSOCK", ps, fixed = TRUE)]
    if (!length(w)) return(c(mem = NA_real_, n = 0))
    rss <- as.numeric(sub("^\\s*([0-9]+).*", "\\1", w))
    c(mem = max(rss, na.rm = TRUE) / 1024, n = length(w))
}

dir_gb <- function(path) {
    if (!dir.exists(path)) return(0)
    out <- suppressWarnings(system2("du", c("-sk", path), stdout = TRUE,
        stderr = FALSE))
    as.numeric(sub("\\s.*", "", out[1L])) / 1024^2
}

write_trace <- function(stage, mem, gb, nw) {
    row <- data.frame(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stage = stage, peak_worker_memory_mb = round(mem, 1),
        cache_size_gb = round(gb, 3), n_workers_seen = nw)
    utils::write.table(row, TRACE, sep = ",", row.names = FALSE,
        col.names = !file.exists(TRACE), append = file.exists(TRACE),
        quote = FALSE)
}

run_stage <- function(stage) {
    cfg <- STAGES[[stage]]
    store <- paste0("_targets_", cfg$tier)
    env <- c(CDA_SIM_TIER = cfg$tier, CDA_SIM_BATCH = BATCH,
        CDA_SIM_WORKERS = WORKERS, OMP_NUM_THREADS = 1,
        OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1,
        VECLIB_MAXIMUM_THREADS = 1)
    cmd <- sprintf(
        'targets::tar_make_future(workers = %d, store = "%s", callr_function = NULL)',
        WORKERS, store)
    logf <- file.path("results", "scratch", paste0("a24-", stage, ".log"))
    dir.create(dirname(logf), recursive = TRUE, showWarnings = FALSE)
    say("stage", stage, "start; store", store)
    t0 <- Sys.time()
    p <- processx::process$new("Rscript", c("-e", cmd), env = c("current", env),
        stdout = logf, stderr = "2>&1", cleanup_tree = TRUE)
    trace_mem <- 0
    trace_nw <- 0
    last_trace <- Sys.time() - 400
    panic <- NULL
    repeat {
        if (!p$is_alive()) break
        p$wait(30000)
        s <- sample_workers()
        gb <- dir_gb(store)
        if (!is.na(s[["mem"]])) trace_mem <- max(trace_mem, s[["mem"]])
        trace_nw <- max(trace_nw, s[["n"]])
        el_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
        tot_h <- as.numeric(difftime(Sys.time(), session_start, units = "hours"))
        if (!is.na(s[["mem"]]) && s[["mem"]] > 3 * 1024) {
            panic <- sprintf("worker memory %.0f MB > 3 GB", s[["mem"]])
        } else if (gb > 10) {
            panic <- sprintf("targets cache %.2f GB > 10 GB", gb)
        } else if (el_h > 1.5 * cfg$budget_h) {
            panic <- sprintf("stage wall %.2f h > 150%% of budget %.2f h",
                el_h, cfg$budget_h)
        } else if (tot_h > SESSION_CAP_H) {
            panic <- sprintf("session wall %.2f h > %d h", tot_h, SESSION_CAP_H)
        }
        if (as.numeric(difftime(Sys.time(), last_trace, units = "secs")) >= 300) {
            write_trace(stage, trace_mem, gb, trace_nw)
            last_trace <- Sys.time()
            trace_mem <- 0
        }
        if (!is.null(panic)) {
            say("PANIC:", panic)
            p$kill_tree()
            break
        }
    }
    write_trace(stage, trace_mem, dir_gb(store), trace_nw)
    wall_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    list(panic = panic, wall_h = wall_h, exit = p$get_exit_status(),
        store = store, log = logf)
}

interrupted <- function(stage, why, chk = NULL, run = NULL) {
    txt <- paste0("# A2-4 interrupted (stub written by the driver)\n\n",
        "Stage: ", stage, "\n\nReason: ", why, "\n\n",
        if (!is.null(chk)) format_status_block(chk, run$wall_h) else "",
        "\nSee `session-logs/session-", DATE, "-a2-4.md` for the status ",
        "blocks of the stages that completed, and ",
        "`a2-4-resource-trace.csv`.\n")
    ## Append: an earlier stub or summary in this file must not be lost.
    cat(txt, file = INTERRUPTED, append = file.exists(INTERRUPTED))
}

if (!file.exists(LOG)) {
    append_log(sprintf("# A2-4 running log (%s)\n\nStatus blocks appended by `simulations/scripts/a24_driver.R` after each stage.\n\n",
        DATE))
}

for (stage in todo) {
    cfg <- STAGES[[stage]]
    if (!is.null(cfg$log)) {
        LOG <- file.path("..", "session-logs", paste0(cfg$log, ".md"))
        INTERRUPTED <- file.path("..", "session-logs",
            paste0(cfg$log, "-interrupted.md"))
    }
    run <- run_stage(stage)
    ## Worker futures are occasionally interrupted (observed once at Stage 1,
    ## 6 of 1,350 branches); targets resumes only the errored branches.
    for (k in 1:3) {
        if (!is.null(run$panic)) break
        nerr0 <- sum(!is.na(tar_meta(store = run$store,
            fields = c("name", "error"))$error))
        if (nerr0 == 0L) break
        say("retry", k, "for", nerr0, "errored targets")
        r2 <- run_stage(stage)
        run$wall_h <- run$wall_h + r2$wall_h
        run$panic <- r2$panic
    }
    if (!is.null(run$panic)) {
        chk <- list(stage = stage, cells_run = NA_integer_, continue = FALSE,
            hard_stop_flags = .flags(), soft_stop_flags = .flags(),
            panic_exit_flags = .add_flag(.flags(), "panic", stage, NA_real_,
                "see detail", run$panic))
        append_log(format_status_block(chk, run$wall_h))
        interrupted(stage, paste("PANIC EXIT:", run$panic), chk, run)
        quit(status = 3L)
    }
    say("stage", stage, "pipeline finished in", sprintf("%.2f h", run$wall_h),
        "exit", run$exit)
    grid <- build_grid(cfg$tier)
    meta <- tar_meta(store = run$store, fields = c("name", "error"))
    nerr <- sum(!is.na(meta$error))
    sim <- tar_read(sim_all, store = run$store)
    ## 12 arm/variant rows per dataset; Sobol surface points run the LMM only
    ## (1 row); sensitivity cells drop the slope-only arm (11 rows).
    expected <- sum(grid$reps * ifelse(grid$family == "power_surface", 1L,
        ifelse(grepl("^sens_", grid$family), 11L, 12L)))
    say("rows", nrow(sim), "expected", expected, "errored targets", nerr)
    if (nrow(sim) != expected || nerr > 0L) {
        why <- sprintf("incomplete: %d rows of %d expected, %d errored targets (see %s)",
            nrow(sim), expected, nerr, run$log)
        append_log(paste0("### Stage ", stage, " INCOMPLETE\n\n", why, "\n\n"))
        interrupted(stage, why, NULL, run)
        quit(status = 4L)
    }
    files <- cfg$writer(sim, grid, cfg$out)
    say("tables written:", length(files))
    res <- a24_attach(sim, grid)
    tr <- utils::read.csv(TRACE)
    tr <- tr[tr$stage == stage, ]
    chk <- check_stage_results(cfg$check %||% stage, res, grid, resource = tr,
        wall_hours = run$wall_h, budget_hours = cfg$budget_h,
        total_wall_hours = as.numeric(difftime(Sys.time(), session_start,
            units = "hours")))
    append_log(format_status_block(chk, run$wall_h,
        extra = paste0("Peak worker memory (MB): ",
            if (nrow(tr)) max(tr$peak_worker_memory_mb, na.rm = TRUE) else NA,
            "; cache size (GB): ",
            if (nrow(tr)) max(tr$cache_size_gb, na.rm = TRUE) else NA,
            "; datasets: ", sum(grid$reps),
            "; files: ", paste(basename(files), collapse = ", "))))
    saveRDS(chk, file.path("results", "scratch", paste0("a24-check-", stage,
        ".rds")))
    rm(sim, res)
    invisible(gc())
    if (!chk$continue) {
        why <- if (nrow(chk$panic_exit_flags)) "PANIC EXIT" else "HARD STOP"
        interrupted(stage, why, chk, run)
        quit(status = 2L)
    }
}
say("all requested stages complete")
