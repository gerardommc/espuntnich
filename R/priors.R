#' @title Generate a list of priors for espuntnich models
#' @description
#' Use a local type of CovMat ppmve or ppmveAsym model to set the prior values for a locallocal model
#' @param model A model object of class ppmve or ppmveAsym
#' @param  ... parameters passed on to priors
#' @return Returns a list with mean and precision for centroids, covariance matrices and intercept
#' @export
#' @method priors


priors <- function(model, ...){
  
  UseMethod(generic = "priors", object = model)  
}