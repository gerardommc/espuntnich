#' @title Learn how to set up priors for a Euclidean distance point process model
#' @description
#' Identify the correct structure to specify the priors for a euclidean distance model
#' @param cent.mean Numeric, with as many elements as covariates in the ppmve model. Must be specified in the same order they appear in "covariate.names".
#' Each value indicates the likely niche centroid coordinates.
#' @param cent.prec Numeric, with as many elements as covariates in the ppmve model. Must be specified in the same order they appear in "covariate.names".
#' Each value indicates the precision (1/variance) of the niche centroid coordinates.
#' @param tau.min Positive numeric, with as many elements as covariates in the ppmve model. Must be specified in the same order they appear in "covariate.names".
#' Each value indicates the minimum value that precision is likely to have.
#' @param tau.max Positive numeric, with as many elements as covariates in the ppmve model. Must be specified in the same order they appear in "covariate.names".
#' Each value indicates the maximum value that precision is likely to have.
#' @param beta.mean A single numeric positive or negative value, inticates the mean of the global intercept.
#' @param beta.prec A single numeric positive, inticates the precision of the global intercept.
#' @examples
#' priors <- euclidPriors(cent.mean = rep(0, 4),
#'                        cent.prec = rep(0.1, 4),
#'                        tau.min = rep(0.01, 4),
#'                        tau.max = rep(100, 4),
#'                        beta.mean = 0,
#'                        beta.prec = 1.0E-4)
#' @export

euclidPiors <- function(cent.mean = NULL,
                        cent.prec = NULL,
                        tau.min = NULL,
                        tau.max = NULL,
                        beta.mean = NULL,
                        beta.prec = NULL){
  
  if(length(cent.mean) != length(cent.prec) | length(tau.min) != length(tau.max) | length(cent.prec) != length(tau.min)){
    stop("Please make sure that cent.mean, cent.prec, tau.min and tau.max contain the same  number of elements")
  }

  if(any(tau.min > tau.max)){
    stop("Please make sure that all values of tau.min are less than the corresponding tau.max values")
  }

  if(any(tau.min < 0) | any(tau.max < 0) | beta.prec < 0){
    stop("Please make sure that all values of tau.min, tau.mx and beta.prec are positive")
  }

  if(length(beta.mean) > 1 | length(beta.prec) > 1){
    stop("Please provide a sigle value for beta.mean and beta.prec")
  }
  
  priors <- list(cent.mean = cent.mean,
                 cent.prec = cent.prec,
                 tau.min = tau.min, 
                 tau.max = tau.max, 
                 beta.mean = beta.mean,
                 beta.prec = beta.prec)

  return(priors)
}