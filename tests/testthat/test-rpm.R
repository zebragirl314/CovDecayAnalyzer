## Tests for rpm.R.

test_that("compute_viral_rpm returns the expected values", {
    rpm <- compute_viral_rpm(
        viral_reads = c(1e4, 5e3, 2e4),
        human_reads = c(1e7, 1e7, 1e6)
    )
    expect_equal(rpm, c(1000, 500, 20000))
})


test_that("compute_viral_rpm handles zero human reads with NA", {
    expect_warning(
        rpm <- compute_viral_rpm(
            viral_reads = c(100, 200),
            human_reads = c(1e6, 0)
        ),
        "Zero human reads"
    )
    expect_equal(rpm[1], 100)
    expect_true(is.na(rpm[2]))
})


test_that("compute_viral_rpm errors on length mismatch", {
    expect_error(
        compute_viral_rpm(viral_reads = 1:3, human_reads = 1:2),
        "same length"
    )
})


test_that("compute_viral_rpm errors on negative values", {
    expect_error(
        compute_viral_rpm(viral_reads = c(-1, 2), human_reads = c(1e6, 1e6)),
        "non-negative"
    )
})


test_that("plot_rpm_by_group returns a ggplot object", {
    set.seed(1)
    rpm_df <- data.frame(
        sample    = paste0("s", 1:12),
        group     = rep(c("WT-GCV", "WT+GCV", "FV-GCV", "FV+GCV"),
                        each = 3),
        viral_rpm = c(rnorm(3, 3000, 500), rnorm(3, 1500, 300),
                      rnorm(3, 2500, 400), rnorm(3, 1200, 250))
    )
    p <- plot_rpm_by_group(rpm_df)
    expect_s3_class(p, "ggplot")
})


test_that("plot_rpm_by_group errors on missing columns", {
    expect_error(
        plot_rpm_by_group(data.frame(x = 1:3, y = 1:3)),
        "required column"
    )
})
