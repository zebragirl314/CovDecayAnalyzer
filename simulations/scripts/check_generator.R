## Ad-hoc generator validation: runs checks 1 to 6 (PLAN.md section 6) and
## prints a report. Run from simulations/:
##   Rscript scripts/check_generator.R
## The testthat version (tests/testthat/test-generator.R) runs the same
## checks for CI-style use. CDA_SIM_CHECK_R raises the calibration and
## positive-control replicate counts.

source("R/setup.R")
load_sim()

results <- list()
for (fn in c("check_parameter_recovery", "check_calibration",
    "check_positive_control", "check_path_equivalence",
    "check_determinism", "check_normalization")) {
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(get(fn)(), error = function(e) {
        list(name = fn, pass = FALSE,
            detail = paste("ERROR:", conditionMessage(e)))
    })
    res$seconds <- proc.time()[["elapsed"]] - t0
    cat(sprintf("[%s] %s (%.0f s)\n    %s\n",
        if (res$pass) "PASS" else "FAIL", res$name, res$seconds, res$detail))
    results[[fn]] <- res
}

n_fail <- sum(!vapply(results, `[[`, logical(1L), "pass"))
cat(sprintf("\n%d of %d checks passed.\n", length(results) - n_fail,
    length(results)))
if (n_fail > 0L) quit(status = 1L)
