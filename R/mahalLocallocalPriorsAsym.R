#' @title Set up priors for an asymmetric Mahalanobis distance point process model
#' @description Identify the correct structure to specify the priors for an asymmetric ellipsoid model
#' @param cent.mean Numeric vector indicating the likely niche centroid coordinates.
#' @param cent.prec Numeric vector indicating the precision of the niche centroid coordinates.
#' @param R A square positive-definite matrix array representing the inverse covariance layout.
#' @param alpha.mean Numeric vector indicating the mean of the asymmetric skewness parameters.
#' @param alpha.prec Numeric vector indicating the precision of the asymmetric skewness parameters.
#' @param beta.mean A single numeric value indicating the mean of the global intercept.
#' @param beta.prec A single numeric positive value indicating the precision of the global intercept.
#' @examples
#' priors <- mahalLocallocalPriorsAsym(cent.mean = rep(0, 4),
#'                                 cent.prec = rep(1.0E-4, 4),
#'                                 R = diag(4),
#'                                 alpha.mean = rep(0, 4),
#'                                 alpha.prec = rep(0.1, 4),
#'                                 beta.mean = 0,
#'                                 beta.prec = 1.0E-4)
#' @export

mahalLocallocalPriorsAsym <- function(cent.mean = NULL,
                                      cent.prec = NULL,
                                      R = NULL,
                                      alpha.mean = NULL,
                                      alpha.prec = NULL,
                                      beta.mean = NULL,
                                      beta.prec = NULL){
  
  if(dim(R)[1] != dim(R)[2] | length(cent.mean) != length(cent.prec) | length(cent.mean) != dim(R)[1]){
    stop("Please make sure that cent.mean, cent.prec, and the dimensions of R match your environmental covariates.")
  }

  if(length(alpha.mean) != length(alpha.prec) | length(alpha.mean) != length(cent.mean)){
    stop("Please make sure that alpha.mean and alpha.prec contain the same number of elements as your covariates.")
  }

  if(any(cent.prec < 0)){
    stop("Please make sure that all values of cent.prec are positive.")
  }

  if(any(alpha.prec < 0)){
    stop("Please make sure that all values of alpha.prec are positive.")
  }

  if(beta.prec < 0){
    stop("Please make sure that beta.prec is positive.")
  }

  if(length(beta.mean) > 1 | length(beta.prec) > 1){
    stop("Please provide a single value for beta.mean and beta.prec.")
  }
  
  priors <- list(cent.mean = cent.mean,
                 cent.prec = cent.prec,
                 R = R,
                 alpha.mean = alpha.mean,
                 alpha.prec = alpha.prec,
                 beta.mean = beta.mean,
                 beta.prec = beta.prec)

  return(priors)
}