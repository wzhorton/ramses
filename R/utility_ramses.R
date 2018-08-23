#### utility_extras.R ####


#' Acceptance Rate Calculator
#'
#' Calculates the acceptance rate of an MCMC chain by looking at the number of repeats.
#'
#' @param chain vector of mcmc values
#'
#' @useDynLib ramses, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @export

acc_rate <- function(chain) {
  n <- length(chain)
  accepts <- sapply(1:(n - 1), function(i) chain[i] != chain[i + 1])
  return(mean(accepts))
}

#' Subsetting Vector By Subject ID Number
#' 
#' Generates a vector used for subsetting given a subject number, 
#' the total number of subjects, and the length of the vector to be subsetted.
#' 
#' @param sub_id and integer specifying what subject is being considered.
#' @param total_sub the total number of subjects accounted for in the vector.
#' @param vec vector to be subsetted.
#' @export

subject_subset <- function(sub_id, total_sub, vec){
  n <- length(vec)
  nper_sub <- n/total_sub
  start <- ((1:total_sub)-1)*nper_sub + 1
  vec[seq(from = start[sub_id], len = nper_sub)]
}
