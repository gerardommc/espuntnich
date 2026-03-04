#' @title Predict method for ppmve models.
#' @description
#' Project a ppmve model onto geographic space.
#' @param object A model object of class ppmve
#' @export
#' @method predict ppmve

predict.ppmve <- function(object, ...){
  UseMethod(generic = "predict.ppmve", object)
}

