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


#' Stack a List of Vectors or Matrices
#'
#' Takes a list of vectors or matrices and returns one vector or matrix along with
#' a vector that specifies the first indices corresponding to the original elements.
#' 
#' @param x list of vectors or matrices with the same number of columns
#' @return a list containing the stacked object and a vector of indices the specify where 
#'   each of the original elements begins.
#' @export

stack <- function(x){
  n <- length(x)
  x <- lapply(x, as.matrix)
  
  lens <- sapply(x, nrow)
  inds <- c(1, sapply(1:n, function(i) 1 + sum(lens[1:i]))[-n])
  out <- list()
  out$stack <- do.call(rbind, x)
  out$inds <- inds
  return(out)
}


#' Unstack a Vector or Matrix
#'
#' Takes a vector or matrice and an index vector and returns a list containing the pieces.
#' 
#' @param x vector or matrix.
#' @param inds vector indicating how to split x. If not provided then it will 
#'   automatically attempt to split into equal parts.
#' @param n indicates the number of elements to split into. 
#'   Only needed when inds is not provided.
#' @return list containing the split elements.
#' @export

unstack <- function(x, inds, n){
  x <- as.matrix(x)
  nr <- nrow(x)
  
  if(missing(inds)){
    if((nper_sub <- nr/n) %% 1 != 0) stop("the length of x is not evenly divided by n")
    inds <- ((1:n)-1)*nper_sub + 1
  }
  
  inds <- c(inds, nr + 1)
  
  return(lapply(1:(length(inds)-1), function(i) x[inds[i]:(inds[i+1] - 1),]))
}

