## Generator validation (PLAN.md section 6). Run with
##   testthat::test_dir("simulations/tests/testthat")
## not devtools::test(). The check logic lives in simulations/R/checks.R and
## is shared with scripts/check_generator.R.

expect_check <- function(res) {
    testthat::expect_true(res$pass, info = res$detail)
    testthat::succeed(res$detail)
}

test_that("1. parameter recovery at large n", {
    expect_check(check_parameter_recovery())
})

test_that("2. idealized calibration: Type I near 0.05, p-values uniform", {
    expect_check(check_calibration())
})

test_that("3. positive control: very large Delta gives power near 1", {
    expect_check(check_positive_control())
})

test_that("4. atom path equals per-base load_all_samples path", {
    expect_check(check_path_equivalence())
})

test_that("5. determinism: replicates regenerate from their seeds", {
    expect_check(check_determinism())
})

test_that("6. normalization is exercised: depth alone changes nothing", {
    expect_check(check_normalization())
})
