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

acc_rate <- function(chain){
  n <- length(chain)
  lag1chain <- chain[-1]
  chain <- chain[-n]
  accepts <- sapply(1:(n-1), function(i) lag1chain[i] != chain[i])
  return(mean(accepts))
}
