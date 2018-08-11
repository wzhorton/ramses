#### compile_ramses.R ####
# This file is only meant to be used for quick recompilation during development.
# For regular installation use devtools::install_github
# Use devtools::install_github(build_vignettes = TRUE) to get html vignettes

library(roxygen2)
library(devtools)

rm(list = ls()) # cleans current environment
gc() # releases memory removed above

file.remove("src/RcppExports.cpp")
file.remove("R/RcppExports.R")


# create NAMESPACE and .Rd files for documentation. Compile C files.
devtools::document("../ramses")

# install preliminary package for function testing and documentation inspection
# setting build_vignettes = FALSE will considerably speed up installation.
# if built with vignettes, use browseVignettes("ramses")
devtools::install("../ramses", build_vignettes = FALSE)
