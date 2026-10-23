## =========================================================================
## Loading the package under test and the simulation code.
##
## The simulation study calls the package's own `bin`/`normalize` code and
## `test_decay_slopes()` rather than reimplementing them, so every entry
## point (scripts, tests, targets workers) needs the package namespace
## loaded. The package is loaded from source with pkgload so the study
## always runs against the working tree, and internals are reachable with
## `:::` exactly as in an installed package.
## =========================================================================

## Locate the simulations/ directory from any working directory inside the
## repository. CDA_SIM_ROOT overrides the search.
sim_root <- function(start = getwd()) {
    env <- Sys.getenv("CDA_SIM_ROOT", unset = "")
    if (nzchar(env)) {
        return(normalizePath(env, mustWork = TRUE))
    }
    dir <- normalizePath(start, mustWork = TRUE)
    repeat {
        if (file.exists(file.path(dir, "_targets.R")) &&
            file.exists(file.path(dir, "PLAN.md"))) {
            return(dir)
        }
        sub <- file.path(dir, "simulations")
        if (file.exists(file.path(sub, "_targets.R"))) {
            return(sub)
        }
        parent <- dirname(dir)
        if (identical(parent, dir)) {
            stop("Could not locate simulations/ from ", start,
                "; set CDA_SIM_ROOT.")
        }
        dir <- parent
    }
}

## Load CovDecayAnalyzer from the repository root (parent of simulations/)
## unless it is already attached or installed in the library.
ensure_package <- function(root = sim_root()) {
    if ("CovDecayAnalyzer" %in% loadedNamespaces()) {
        return(invisible(TRUE))
    }
    pkg_dir <- dirname(root)
    pkgload::load_all(pkg_dir, quiet = TRUE, helpers = FALSE)
    invisible(TRUE)
}

## Source every simulation script in R/ (definitions only; no side effects).
load_sim <- function(root = sim_root()) {
    ensure_package(root)
    files <- sort(list.files(file.path(root, "R"), pattern = "\\.R$",
        full.names = TRUE))
    files <- files[basename(files) != "setup.R"]
    for (f in files) {
        sys.source(f, envir = globalenv())
    }
    invisible(files)
}
