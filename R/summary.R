#' @title Summary methods for ppmve objects
#' @description 
#' Prints the estimated coefficients of a ppmve model
#' @param model A model object of class ppmve
#' @return A summary of the posterior samples for the fitted model
#' @export

summary <- function(model){
  UseMethod(generic = "summary", model)
}