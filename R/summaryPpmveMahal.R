#' @title Print summary statistics for a ppmve mahalanobis distance model.
#' @description
#' This function prints the estimated coefficients of a ppmve model and calculates
#' probability values for each of the parameter posteriors, relevant
#' to the type of analysis for which the method has been designed
#' @param model = A model object of class ppmve
#' @return A summary of the posterior samples for the fitted model
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
#' 
#' # Model with random sampling
#' m <- ppmve(points = pr,
#'            covariates = r,
#'            covariate.names = names(r),
#'            CovMat = "local",
#'            Distance = "mahalanobis",
#'            no.bkgd = 5000,
#'            niter = 10000,
#'            nthin = 9,
#'            nburnin = 1000,
#'            chains = 1)
#' 
#' summary(m)
#' @export
#' @method summary.ppmve mahalanobis

summary.ppmve.mahalanobis <- function(model){
  summary(model$model$samples)
}


