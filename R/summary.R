#' @title Print summary statistics for a ppmve model.
#' @description
#' This function prints the estimated coefficients of a ppmve model and calculates
#' probability values for each of the parameter posteriors, relevant
#' to the type of analysis for which the method has been designed
#' @param model = A model object of class ppmve
#' @return A summary of the posterior samples for the fitted model


summary.ppmve <- function(model = NULL){
  UseMethod("summary.ppmve", model)
}

summary.ppmve.mahalanobis <- function(model = NULL){
  summary(model$model$samples)
}

summary.ppmve.euclidean <- function(model = NULL){
  summary(model$model$samples)
}

