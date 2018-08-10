#### compile_ramses.R ####
# This file is only meant to be used for quick recompilation during development.
# For regular installation use devtools::install_github

library(roxygen2)
library(devtools)

rm(list = ls()) # cleans current environment
gc() # releases memory removed above

file.remove("src/RcppExports.cpp")
file.remove("R/RcppExports.R")


# create NAMESPACE and .Rd files for documentation. Compile C files.
devtools::document("../ramses")

# install preliminary package for function testing and documentation inspection
devtools::install("../ramses")
