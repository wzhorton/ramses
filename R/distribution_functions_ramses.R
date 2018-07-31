#### distribution_functions_ramses.R ####

#' Inverse Gamma Density Function
#' 
#' Computes the (log) density of the inverse gamma distribution using 
#' either the scale or rate parametrization.
#' 
#' @param x vector of positive values.
#' @param shape,scale shape and scale parameters. Must be stricly positive.
#' @param rate alternative way to specify scale.
#' @param log logical; if TRUE, density calculations are computed on the log scale.
#' @export

dinvgamma <- function(x, shape, rate, scale = 1/rate, log = FALSE){
  if(shape <= 0 || scale <= 0 || x <= 0) stop("Shape, rate, scale, and x 
                                              must all be positive")
  a <- shape
  b <- scale
  out <- -a*log(b)-lgamma(a)+(-a-1)*log(x)-1/(x*b)
  if(log == FALSE) out <- exp(out)
  return(out)
}

#' Random Multivariate Normal Generator
#' 
#' Generates a normally distributed vector given a mean vector and 
#' either a covariance matrix or a precision matrix. Does not support 
#' generating multiple vectors at once.
#' 
#' @param mu mean vector.
#' @param cov covariance matrix.
#' @param prec precision matrix.
#' 
#' @return A randomly generated vector with the same length as mu.
#' @export

rmnorm <- function(mu, cov, prec){
  if(missing(mu)) stop("Provide a mean vector")
  if(!xor(missing(prec), missing(cov))){ 
    stop("Provide either Precision or Covariance, but not both")
  }
  if(missing(prec)) out <- as.numeric(mu + t(chol(cov))%*%rnorm(length(mu)))
  else if(missing(cov)) out <- mu +  backsolve(chol(prec), rnorm(length(mu)))
  return(out)
}

#' Multivariate Normal Density Function
#' 
#' Computes the (log) density of the multivariate normal distribution 
#' using either the covariance or precision parametrization.
#' 
#' @param y vector of values.
#' @param mu mean vector.
#' @param cov covariance matrix.
#' @param prec precision matrix.
#' @param log logical; if TRUE, density calculations are computed on the log scale.
#' @export

dmnorm <- function(y, mu, cov, prec, log = FALSE){
  if(!xor(missing(cov),missing(prec))) stop("Provide either cov or prec, but not both")
  n <- length(y)
  
  if(missing(prec)){
    out <- as.numeric(-n/2*log(2*pi)-.5*log(det(as.matrix(cov)))-
                        .5*t(y - mu)%*%chol2inv(chol(cov))%*%(y-mu))
  }
  else if(missing(cov)){
    out <- as.numeric(-n/2*log(2*pi)+.5*log(det(as.matrix(prec)))-
                        .5*t(y - mu)%*%prec%*%(y-mu))
  }
  
  if(log == FALSE) out <- exp(out)
  return(out)
}




