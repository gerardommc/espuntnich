#' @title Predict method for ppmveAsym models.
#' @description
#' Project an asymmetric ppmve model onto geographic space.
#' @param object A model object of class ppmveAsym
#' @param newdata A SpatRaster used for projecting the ppmveAsym model
#' @param probs A numeric value indicating the percentiles to be shown in the predictions
#' @export
#' @method predict ppmveAsym

predict.ppmveAsym <- function(object, newdata, probs){
  UseMethod(generic = "ppmveAsym", object)
}