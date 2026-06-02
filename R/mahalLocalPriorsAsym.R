#' @title Set up priors for a Mahalanobis distance point process model with local background covariance tracking (Asymmetric)
#' @description Identify the correct structure to specify the priors for an asymmetric ellipsoid model using the local background covariance approach.
#' @param cent.mean Numeric vector with as many elements as covariates in the model, indicating the likely niche centroid coordinates.
#' @param cent.prec Numeric vector with as many elements as covariates in the model, indicating the precision of the niche centroid coordinates.
#' @param mu.b.mean Numeric vector with as many elements as covariates in the model, indicating the mean of the background environmental distribution.
#' @param mu.b.prec Numeric vector with as many elements as covariates in the model, indicating the precision of the background environmental distribution.
#' @param R A square positive-definite matrix representing the scale matrix for the Wishart prior on the precision matrix.
#' @param alpha.mean Numeric vector with as many elements as covariates in the model, indicating the mean of the asymmetric skewness parameters ($\alpha$).
#' @param alpha.prec Numeric vector with as many elements as covariates in the model, indicating the precision of the asymmetric skewness parameters ($\alpha$).
#' @param beta.mean A single numeric value indicating the mean of the global intercept.
#' @param beta.prec A single numeric positive value indicating the precision of the global intercept.
#' @examples
#' priors <- mahalLocalPriorsAsym(cent.mean = rep(0, 4),
#'                                cent.prec = rep(1.0E-4, 4),
#'                                mu.b.mean = rep(0, 4),
#'                                mu.b.prec = rep(1.0E-4, 4),
#'                                R = diag(4),
#'                                alpha.mean = rep(0, 4),
#'                                alpha.prec = rep(0.1, 4),
#'                                beta.mean = 0,
#'                                beta.prec = 1.0E-4)
#' @export

mahalLocalPriorsAsym <- function(cent.mean = NULL,
                                 cent.prec = NULL,
                                 mu.b.mean = NULL,
                                 mu.b.prec = NULL,
                                 R = NULL,
                                 alpha.mean = NULL,
                                 alpha.prec = NULL,
                                 beta.mean = NULL,
                                 beta.prec = NULL){
  
  # Check dimension compatibility for multivariate parameters
  n_clim <- dim(R)[1]
  if(dim(R)[1] != dim(R)[2]){
    stop("The scale matrix R must be a square matrix.")
  }
  
  if(length(cent.mean) != n_clim | length(cent.prec) != n_clim){
    stop("Please ensure that cent.mean and cent.prec have lengths equal to the dimensions of matrix R.")
  }
  
  if(length(mu.b.mean) != n_clim | length(mu.b.prec) != n_clim){
    stop("Please ensure that mu.b.mean and mu.b.prec have lengths equal to the dimensions of matrix R.")
  }
  
  if(length(alpha.mean) != n_clim | length(alpha.prec) != n_clim){
    stop("Please ensure that alpha.mean and alpha.prec have lengths equal to the dimensions of matrix R.")
  }

  # Precision positivity validations
  if(any(cent.prec < 0)){
    stop("Please make sure that all values of cent.prec are positive.")
  }
  
  if(any(mu.b.prec < 0)){
    stop("Please make sure that all values of mu.b.prec are positive.")
  }

  if(any(alpha.prec < 0)){
    stop("Please make sure that all values of alpha.prec are positive.")
  }

  if(beta.prec < 0){
    stop("Please make sure that beta.prec is positive.")
  }

  # Intercept length validations
  if(length(beta.mean) > 1 | length(beta.prec) > 1){
    stop("Please provide a single value for beta.mean and beta.prec.")
  }
  
  # Package priors list
  priors <- list(cent.mean  = cent.mean,
                 cent.prec  = cent.prec,
                 mu.b.mean  = mu.b.mean,
                 mu.b.prec  = mu.b.prec,
                 R          = R,
                 alpha.mean = alpha.mean,
                 alpha.prec = alpha.prec,
                 beta.mean  = beta.mean,
                 beta.prec  = beta.prec)

  return(priors)
}