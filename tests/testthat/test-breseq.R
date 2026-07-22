## Tests for breseq.R.

test_that("breseq_coverage_paths builds default paths correctly", {
    paths <- breseq_coverage_paths(
        sample_names = c("A1", "B2"),
        base_dir     = "/data/out"
    )
    expect_equal(paths, c(
        "/data/out/A1/08_mutation_identification/Exported.coverage.tab",
        "/data/out/B2/08_mutation_identification/Exported.coverage.tab"
    ))
})


test_that("breseq_coverage_paths honors custom subdir and filename", {
    paths <- breseq_coverage_paths(
        sample_names = "S1",
        base_dir     = "/out",
        subdir       = "07_error_calibration",
        filename     = "coverage.txt"
    )
    expect_equal(paths, "/out/S1/07_error_calibration/coverage.txt")
})


test_that("breseq_coverage_paths preserves sample order", {
    ss <- c("z", "a", "m", "b")
    paths <- breseq_coverage_paths(ss, base_dir = "/x")
    expect_equal(basename(dirname(dirname(paths))), ss)
})


test_that("breseq_coverage_paths rejects empty or NA sample vectors", {
    expect_error(
        breseq_coverage_paths(character(0), base_dir = "/x"),
        "non-empty"
    )
    expect_error(
        breseq_coverage_paths(c("A", NA_character_), base_dir = "/x"),
        "NA"
    )
})
