#' @title Print summary statistics for a ppmve euclidean distance model.
#' @description
#' This function prints the estimated coefficients of a ppmve model and calculates
#' probability values for each of the parameter posteriors, relevant
#' to the type of analysis for which the method has been designed
#' @param model A model object of class ppmve
#' @return A summary of the posterior samples for the fitted model
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espuntnich") |> terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espuntnich") |> read.csv()
#' 
#' m <- ppmve(points = p,
#'            covariates = r,
#'            covariate.names = names(r),
#'            Distance = "euclidean",
#'            no.bkgd = 5000,
#'            niter = 10000,
#'            nthin = 9,
#'            nburnin = 1000,
#'            chains = 1)
#' 
#' summary.ppmve(m)
#' @export
#' @method summary.ppmve euclidean

summary.ppmve.euclidean <- function(model){
  summary(model$model$samples)
}


