#### utility_general.R ####

#' Spline Curve Interpolation
#'
#' Interpolates out equally spaced values by fitting a cubic spline on given points and
#' matching the endpoints. Often used to redefine a curve with warped time axis on an
#' even time basis. Note that quality of interpolation decreases with inadequate density
#' of defining points.
#'
#' @param x,y coordinate vectors on which to fit a cubic spline.
#' @param nout number of points to interpolate out
#' @export

interp_spline <- function(x, y, nout = length(y)) {
  ind_out <- seq(min(x), max(x), len = nout)
  spfit <- splinefun(x, y)
  return(spfit(ind_out))
}

#' Second Order Penalty Matrix
#'
#' Creates a second order penalty matrix used in fitting penalized b-splines.
#' Courtesy of CITATION NEEDED. Note that this function produces singular matrices.
#'
#' @param dim dimension of output matrix. Must be at least four.
#' @export

K2 <- function(dim) {
  tmp <- list(c(-2, rep(-4, dim - 3), -2), rep(1, dim))
  K <- Matrix::bandSparse(dim, k = -c(1:2), diag = tmp, symmetric = TRUE)
  diag(K) <- c(1, 5, rep(6, dim - 4), 5, 1)
  K <- matrix(as.numeric(K), nrow = dim, byrow = TRUE)
  K
}
