#' @title Print summary statistics for a ppmve model.
#' @description
#' This function prints the estimated coefficients of a ppmve model and calculates
#' probability values for each of the parameter posteriors, relevant
#' to the type of analysis for which the method has been designed
#' @param model A model object of class ppmve
#' @return A summary of the posterior samples for the fitted model
#' @export
#' @method summary ppmve

summary.ppmve <- function(model){
  UseMethod(generic = "summary.ppmve", model)
}


