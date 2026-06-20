#' @title Predict method for ppmve models.
#' @description
#' Project a ppmve model onto geographic space.
#' @param object A model object of class ppmve
#' @param newdata A SpatRaster used for projecting the ppmve model
#' @param probs A numeric value indicating the percentiles to be shown in the predictions
#' @export
#' @method predict ppmve

predict.ppmve <- function(object, newdata, probs){
  UseMethod(generic = "ppmve", object)
}

