## =========================================================================
## Internal input-validation helpers. Not exported.
## =========================================================================

#' @keywords internal
#' @noRd
.validate_string <- function(x, name) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        stop(sprintf("`%s` must be a single non-empty character string.", name))
    }
    invisible(TRUE)
}

#' @keywords internal
#' @noRd
.validate_positive_integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x <= 0 || x != as.integer(x)) {
        stop(sprintf("`%s` must be a single positive integer.", name))
    }
    invisible(TRUE)
}

#' @keywords internal
#' @noRd
.validate_position <- function(x, name, genome_length = NULL) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
        stop(sprintf(
            "`%s` must be a single non-negative finite number.", name
        ))
    }
    if (!is.null(genome_length) && x > genome_length) {
        stop(sprintf(
            "`%s` (%g) exceeds `genome_length` (%g).",
            name, x, genome_length
        ))
    }
    invisible(TRUE)
}
