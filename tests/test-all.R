## This script is only run during R CMD check, so we can set an environment
## variable that will only run tests during R CMD check (or `devtools::check()`)
## and not during `devtools::test()`.

## Thus in the tests, we can request the NEXUS files that are stored in the
## `inst/` folder, but during the checks, we test the files that have been
## installed (using the `system.file()` function).

library(testthat)

Sys.setenv("R_CMD_CHECK" = "true")

test_check("phylobase")
