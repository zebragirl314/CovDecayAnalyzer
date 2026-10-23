## Loads the package under test and the simulation code, whatever the
## working directory of the test run (repo root or simulations/).
.sim_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
    mustWork = TRUE)
source(file.path(.sim_root, "R", "setup.R"))
Sys.setenv(CDA_SIM_ROOT = .sim_root)
load_sim(.sim_root)
