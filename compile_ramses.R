#### compile_ramses.R ####
#only useful during development

library(roxygen2)
library(devtools)

#create NAMESPACE and .Rd files for documentation
document("../ramses")

#install preliminary package for function testing and documentation inspection
install("../ramses")

#notes for next time. Some of the documentation titles are messed up. No armadillo uses yet. 
# try to look for some in the current functions. Also look into compilation needs.
# manual did not compile into pdf. Need vignnets.
#newlines within description for update models.