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
  accepts <- sapply(1:(n-1), function(i) chain[i] != chain[i+1])
  return(mean(accepts))
}
