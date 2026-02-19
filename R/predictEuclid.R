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
#'            Distance = "euclidean",
#'            no.bkgd = 5000,
#'            niter = 10000,
#'            thin = 9,
#'            nburnin = 1000,
#'            chains = 1) 
#' 
#' predictions <- predict(object = m, nedata = r, probs = c(0.0275, 0.5, 0.975))
#' 
#' plot(predictions)
#' }
#' @export

predict.ppmve.euclidean <- function(object = NULL, newdata = NULL, probs = 0.5){

  if(class(newdata) != "SpatRaster"){
    stop("Please provide a SpatRaster object")
  }

  if(class(newdata) == "SpatRaster"){
    cov.df <- newdata |> as.data.frame(xy = T)
    cov.df1 <- cov.df |> subset(select = c(object$call$covariates)) |> as.matrix()
  }

  if(object$call$chains == 1){
    mcl <- object$model$samples |> coda::as.mcmc()
  }

  if(object$call$chains > 1){
    mcl <- object$model$samples |> coda::as.mcmc.list()
    mc <- mcl[[1]]
    for(i in 2:object$call$chains){
      mc <- rbind(mc, mcl[[i]])
    }

    mcl <- mc |> coda::as.mcmc()
  }

  if(length(probs) > 1){
    coefs <- lapply(probs, function(x){coda::HPDinterval(mcl, prob = x)})
  }

  if(length(probs) == 1){
    coefs <- coda::HPDinterval(mcl, prob = probs)
  }

  if(length(probs) > 1){

      p <- foreach::foreach(i = seq_along(probs), .combine = cbind) %do% {

        mu <- coefs[[i]][paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
        beta <- coefs[[i]]["beta", ] |> mean()
        tau.mat <- coefs[[i]][paste0("tau.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans() |> diag()

        md <- mahalanobis(cov.df1, center = mu, cov = tau.mat)
        md.ex <- exp(beta-md/2)

        return(md.ex)

      }

        preds <- data.frame(cov.df[, c("x", "y")], p) |> terra::rast()
        names(preds) <- paste0("Prob.", probs)
        crs(preds) <- terra::crs(newdata)
  }


  if(length(probs) == 1){

      mu <- coefs[paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
      beta <- coefs["beta", ] |> mean()
      tau.mat <- coefs[paste0("tau.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans() |> diag()

      md <- mahalanobis(cov.df1, center = mu, cov = tau.mat)
      md.ex <- exp(beta-md/2)

      preds <- data.frame(cov.df[, c("x", "y")], md.ex) |> terra::rast()
      names(preds) <- paste0("Prob.", probs)
      terra::crs(preds) <- terra::crs(newdata)
  }

  return(preds)
}

