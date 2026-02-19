#' @title Print summary statistics for a ppmve model.
#' @description
#' This function prints the estimated coefficients of a ppmve model and calculates
#' probability values for each of the parameter posteriors, relevant
#' to the type of analysis for which the method has been designed
#' @param model = A model object of class ppmve
#' @return A summary of the posterior samples for the fitted model
#' @examples
#' \dontrun{
#' r <- terra::rast(paste0("bio",c(1, 2, 12, 17), ".tif")) |> scale()
#' 
#' p <- read.csv("points.csv")
#' 
#' m <- ppmve(points = p,
#'            covariates = r,
#'            covariate.names = names(r),
#'            CovMat = "local",
#'            Distance = "mahalanobis",
#'            no.bkgd = 5000,
#'            niter = 10000,
#'            thin = 9,
#'            nburnin = 1000,
#'            chains = 1)
#' 
#' summary(model$model$samples)
#' }
#' @export

summary.ppmve <- function(model = NULL){
  summary(model = model$model$samples)
}


