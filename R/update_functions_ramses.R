#### update_functions_ramses.R ####

#' Conjugate Normal - Inverse Gamma Update
#'
#' Generates a value from the posterior distribution in the case where there
#' is a multivariate normal likelihood and an inverse gamma prior.\cr \cr
#' Argument model:\cr
#' y ~ Nn(mu, sig2*R)\cr
#' sig2 ~ IG(a,b)\cr
#' @param y vector of values at the likelihood level.
#' @param a prior shape value for inverse gamma.
#' @param b prior SCALE value for inverse gamma.
#' @param mu mean vector for multivariate normal.
#' @param R correlation matrix for multivariate normal.
#' @param R_inv scaled precision matrix, alternative specification to R.
#' @export

update_normal_invgamma <- function(y, a, b, mu, R, R_inv) {
  if (!xor(missing(R), missing(R_inv))) stop("Provide either R or R_inv, but not both")
  if (missing(R_inv)) R_inv <- chol2inv(chol(R))
  return(1 / rgamma(1, .5 * length(y) + a, 
                    rate = .5 * t(y - mu) %*% R_inv %*% (y - mu) + 1 / b))
}

#' Conjugate Multivariate Normal - Multivariate Normal Update
#'
#' Generates a value from the posterior distribution in the case where
#' there is a multivariate normal likelihood and a multivariate normal prior.\cr \cr
#' Argument Model:\cr
#' y ~ Nn(X*beta, Sig)\cr
#' beta ~ Np(mu, V)\cr
#' @param y vector of values at the likelihood level.
#' @param X fixed design matrix in likelihood.
#' @param mu prior mean vector.
#' @param Sig,Sig_inv likelihood covariance/precision matrix.
#' @param V,V_inv prior covariance/precision matrix.
#' @export

update_normal_normal <- function(y, X, mu, Sig, V, Sig_inv, V_inv) {
  if (!xor(missing(Sig), missing(Sig_inv))) {
    stop("Provide either Sig or Sig_inv, but not both")
  }
  if (!xor(missing(V), missing(V_inv))) stop("Provide either V or V_inv, but not both")
  if (missing(Sig_inv)) Sig_inv <- chol2inv(chol(Sig))
  if (missing(V_inv)) V_inv <- chol2inv(chol(V))

  vv <- t(X) %*% Sig_inv %*% X + V_inv
  vterm <- chol2inv(chol(vv))
  return(rmnorm(vterm %*% (t(X) %*% (Sig_inv %*% y) + V_inv %*% mu), prec = vv))
}

#' Gaussian Process Update
#'
#' Updates and returns the mean vector for a gaussian process given evaluation points,
#' observed data, a mean function, and a covariance function.\cr \cr
#' Argument model:\cr
#' y(time) ~ GP(mnfun(x),covfun(x1-x2))
#'
#' @param x,y coordinate vectors for the observed data points.
#' @param time locations to evaluate the curve at.
#' @param mnfun mean function that takes a single argument.
#'   If using a constant like 0 use function(z)\{0\}.
#' @param covfun covariance function that takes a single argument that represents
#'   an absolute distance. Distances are internally calculated using fields::rdist.
#' @param random logical; FALSE will return the updated mean function. 
#'   TRUE will generate a random curve.
#' @return a vector corresponding to the time variable that represents the updated
#'   mean vector.
#' @export

update_gaussian_process <- function(x, y, time, mnfun, covfun, random = FALSE) {
  R12 <- covfun(fields::rdist(time, x))
  R22 <- covfun(fields::rdist(x))
  mu1 <- mnfun(time)
  mu2 <- mnfun(x)
  up_mean <- .update_gp_mean_C(y, R12, R22, mu1, mu2)
  if(random == FALSE){
    out <- up_mean
  } else {
    R11 <- covfun(fields::rdist(time))
    up_var <- .update_gp_var_C(R11, R12, R22)
    out <- rmnorm(up_mean, cov = up_var)
  }
  return(out)
}

#' Update Step Routine
#'
#' Excecutes one full cycle of an MCMC iteration within the RAMSES sampling
#' structure, updating each parameter in the parm_state list with the update
#' function found in the updates list.
#' @param parm_state list that names parameters and contains current values.
#' @param fixed list of named values that remain fixed in the model.
#' @param updates list of update functions with names matching those in parm_state.
#' @param parm_names vector naming the parameters as the appear in order in both
#'   updates and parm_state.
#' @return A list in the form of parm_state, but with updated values. To be used
#'   within the RAMSES sampler function.
#' @export

update_step <- function(parm_state, fixed, updates, parm_names = names(updates)) {
  for (p in parm_names) {
    parm_state[[p]] <- updates[[p]](parm_state, fixed)
  }
  return(parm_state)
}

#' RAMSES Sampling Structure
#'
#' Function that fits a heirarchical model using MCMC. It requires a list of initial
#' parameter values, a list of needed fixed values, and a list of update functions,
#' all of which are tied together by matching element names.
#' @param inits list of named parameter initial values.
#' @param fixed list of named values that remain fixed in the model.
#' @param updates list of update functions with names matching those in parm_state.
#' @param niter,nburn number of MCMC iterates and burnin.
#' @param parms_to_save character vector indicating which parameters to save.
#'         Defaults to all. Use this to allow for auxiliary onjects like matrices.
#' @return A matrix where each row is an MCMC iterate and each column is a
#'   parameter or part of a parameter vector.
#' @export

sampler <- function(inits, fixed, updates, niter, nburn, 
                    parms_to_save = names(inits), progress = FALSE) {
  
  cat("Progress:  0 %")
  pnames <- names(unlist(inits[parms_to_save]))
  parameter_state_list <- inits
  nrun <- niter + nburn
  save_matrix <- matrix(NA,
    ncol = length(pnames), nrow = nrun,
    dimnames = list(NULL, pnames)
  )
  save_matrix[1, ] <- unlist(inits[parms_to_save])

  for (i in 2:nrun) {
    out <- update_step(parameter_state_list, fixed, updates)
    parameter_state_list <- out
    save_matrix[i, ] <- unlist(out[parms_to_save])
    rm(out)
    if(progress == TRUE){
      if(((i %% round(nrun*.01)) == 0) && i/nrun < 0.1) cat("\b\b\b\b",round(i/nrun*100),"%")
      if(((i %% round(nrun*.01)) == 0) && i/nrun >= 0.1) cat("\b\b\b\b\b",round(i/nrun*100),"%")
    }
  }
  cat("\nSampling Complete\n")
  return(save_matrix[-c(1:nburn), ])
}

#' Metropolis-Hastings Symmetric Update
#'
#' Function that performs a M-H update for a named parameter given the RAMSES
#' structural elements and loglikelihood/logprior functions. Candidate values
#' are drawn from a normal distribution, but it does allow for upper and lower bounds.
#' @param current_value current value of the parameter.
#' @param candidate_value optional; specifies the candidate or proposed value. If 
#'  left out, a random normal vector is generated as the proposal.
#' @param loglik_func,logprior_func log likelihood/prior function that takes
#'  the current value, parm_state, and fixed. It is possible to let
#'  the current value be simply extracted from parm_state.
#' @param parm_state list that names parameters and contains current values.
#' @param fixed list of named values that remain fixed in the model.
#' @param tune tuning parameter (standard deviation) in the normal proposal distribution.
#' @param lower,upper bounds for the parameter to restrict proposed values for 
#'  univariate parameters. Invalid values are reflected back into the domain. 
#'  If both bounds are provided, make sure the tune is small so that reflections won't 
#'  also be invalid.
#' @export

update_metropolis <- function(current_value, candidate_value, loglik_func, logprior_func, parm_state, fixed,
                              tune = .01, lower = NULL, upper = NULL) {
  if(missing(candidate_value)){
    candidate_value <- rmnorm(current_value, cov = tune * diag(length(current_value)))
  }

  if (!is.null(lower) && candidate_value <= lower) {
    candidate_value <- abs(candidate_value - lower) + lower
  }
  if (!is.null(upper) && candidate_value >= upper) {
    candidate_value <- upper - abs(candidate_value - upper)
  }
  if (!is.null(lower) && !is.null(upper) &&
    abs(candidate_value - lower) > (upper - lower) &&
    abs(candidate_value - upper) > (upper - lower)) {
    out <- current_value
  }
  else {
    logbottom <- loglik_func(current_value, parm_state, fixed) + 
      logprior_func(current_value, parm_state, fixed)
    
    logtop <- loglik_func(candidate_value, parm_state, fixed) + 
      logprior_func(candidate_value, parm_state, fixed)
    if (log(runif(1)) < logtop - logbottom) {
      out <- candidate_value
    } else {
      out <- current_value
    }
  }
  return(out)
}
