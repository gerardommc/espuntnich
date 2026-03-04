#' @title Summary method for a ppmve model posterior coefficients.
#' @description
#' Print the estimated coefficients of a ppmve model and calculates
#' probability values for each of the parameter posteriors.
#' @param model A model object of class ppmve
#' @return A summary of the posterior samples for the fitted model
#' @export
#' @method summary ppmve

summary.ppmve <- function(model){
  UseMethod(generic = "summary.ppmve", model)
}


