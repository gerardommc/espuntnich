#' @title Summary method for a ppmveAsym model posterior coefficients.
#' @description
#' Print the estimated coefficients of an asymmetric ppmve model and calculates
#' probability values for each of the parameter posteriors.
#' @param model A model object of class ppmveAsym
#' @return A summary of the posterior samples for the fitted model
#' @export
#' @method summary ppmveAsym

summary.ppmveAsym <- function(model){
  UseMethod(generic = "summary.ppmveAsym", model)
}