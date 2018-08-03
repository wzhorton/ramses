#### compile_ramses.R ####
#only useful during development, to install use devtools::install_github

library(roxygen2)
library(devtools)

#create NAMESPACE and .Rd files for documentation
devtools::document("../ramses")

#install preliminary package for function testing and documentation inspection
devtools::install("../ramses")