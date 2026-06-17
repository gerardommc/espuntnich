#' @title Print summary statistics for an asymmetric ppmve mahalanobis distance model.
#' @description
#' This function prints the estimated coefficients of an asymmetric ppmve model and calculates
#' probability values for each of the parameter posteriors, relevant
#' to the type of analysis for which the method has been designed.
#' @param model A model object of class ppmveAsym
#' @return A summary of the posterior samples for the fitted asymmetric model
#' @examples
#' 
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espuntnich") |> terra::rast() |> scale()
#' p <- system.file("extdata", "pointsAsym.csv", package = "espuntnich") |> read.csv()
#' 
#' m <- ppmveAsym(points = p,
#'                covariates = r,
#'                covariate.names = names(r),
#'                CovMat = "local",
#'                Distance = "mahalanobis",
#'                no.bkgd = 5000,
#'                niter = 10000,
#'                nthin = 9,
#'                nburnin = 1000,
#'                chains = 1)
#' 
#' summary(m)
#' 
#' @export
#' @method summary.ppmveAsym mahalanobis

summary.ppmveAsym.mahalanobis <- function(model){
  summary(model$model$samples)
}