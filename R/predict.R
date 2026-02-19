#' @title Generate predictions from a ppmve model.
#' @description
#' Project a ppmve model onto geographic space. To do this, the user should profide
#' a fitted model, and the covariates in SpatRaster format with the same name as
#' used to fit the model.
#' @param model = A model object of class ppmve
#' @param newdata = A SpatRaster object with covariate names contained in the covariate.names argument of ppmve function
#' @param probs = The posterior probability quantiles to be returned by predict.ppmve
#' @return Returns a single or multiple band SpatRaster object, representing point intensity as a function of distance to the estimated centroids
#' @examples
#' \dontrun{
#' r <- terra::rast(paste0("bio",c(1, 2, 12, 17), ".tif")) |> scale()
#' 
#' p <- read.csv("points.csv")
#' 
#' m <- ppmve(points = p,
#'            covariates = r,
#'            covariate.names = names(r),
#'            CovMat = "local",
#'            Distance = "mahalanobis",
#'            no.bkgd = 5000,
#'            niter = 10000,
#'            thin = 9,
#'            nburnin = 1000,
#'            chains = 1)
#' }
#' 
#' predictions <- predict(model = m, nedata = r, probs = c(0.0275, 0.5, 0.975))
#' 
#' plot(predictions)
#' @export

predict.ppmve <- function(model = NULL, newdata = NULL, probs = 0.5){
  UseMethod("predict.ppmve", model)
}
