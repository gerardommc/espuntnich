#' @title Generate predictions from a ppmve model.
#' @description
#' Project a ppmve model onto geographic space. To do this, the user should profide
#' a fitted model, and the covariates in SpatRaster format with the same name as
#' used to fit the model.
#' @param model = A model object of class ppmve
#' @param newdata = A SpatRaster object with covariate names contained in the covariate.names argument of ppmve function
#' @param probs = The posterior probability quantiles to be returned by predict.ppmve
#' @return Returns a single or multiple band SpatRaster object, representing point intensity as a function of distance to the estimated centroids
#' @export

predict.ppmve <- function(object, ...){
  NextMethod("predict.ppmve", object)
}
