#' @title Generate predictions from an asymmetric ppmve model.
#' @description
#' Project an asymmetric ppmve model onto geographic space. To do this, the user should provide
#' a fitted asymmetric model, and the covariates in SpatRaster format with the same name as
#' used to fit the model.
#' @param object A model object of class ppmveAsym
#' @param newdata A SpatRaster object with covariate names contained in the covariate.names argument of ppmveAsym function
#' @param probs The posterior probability quantiles to be returned by predict.ppmveAsym
#' @return Returns a single or multiple band SpatRaster object, representing point intensity as a function of asymmetric distance to the estimated centroids
#' @export
#' @method predict.ppmveAsym mahalanobis

predict.ppmveAsym.mahalanobis <- function(object,
                                          newdata = NULL,
                                          probs = c(0.0275, 0.5, 0.975)){
  
  `%do%` <- foreach::`%do%`

  if(is.null(newdata)){
    stop("Please provide a valid data.frame or SpatRaster object")
  }

  if(!inherits(newdata, "SpatRaster")){
    stop("Please provide a SpatRaster object")
  }
  
  if(inherits(newdata, "SpatRaster")){
    cov.df <- newdata |> as.data.frame(xy = TRUE)
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

  tau.names <- expand.grid(i = 1:ncol(cov.df1), j = 1:ncol(cov.df1))

  if(length(probs) > 1){

    p <- foreach::foreach(i = seq_along(probs), .combine = cbind) %do% {

      mu <- coefs[[i]][paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
      beta <- coefs[[i]]["beta", ] |> mean()
      tau <- coefs[[i]][paste0("tau.pres[", tau.names$i, ", ", tau.names$j, "]"), ] |> rowMeans()
      tau.mat <- matrix(tau, nrow = ncol(cov.df1), ncol = ncol(cov.df1))

      # 1. Symmetric distance component
      md <- stats::mahalanobis(cov.df1, center = mu, cov = tau.mat)
      
      # 2. Asymmetric skewness component extraction & evaluation
      alpha <- coefs[[i]][paste0("alpha[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
      X_cent <- sweep(cov.df1, 2, mu, "-")
      skew.val <- X_cent %*% alpha
      
      # 3. Combine using the Skew-Normal probability multiplier
      md.ex <- exp(beta - md/2) * 2 * stats::pnorm(skew.val)

      return(md.ex)
    }

    preds <- data.frame(cov.df[, c("x", "y")], p) |> terra::rast()
    names(preds) <- paste0("Prob.", probs)
    terra::crs(preds) <- terra::crs(newdata)
  }

  if(length(probs) == 1){

    mu <- coefs[paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
    beta <- coefs["beta", ] |> mean()
    tau <- coefs[paste0("tau.pres[", tau.names$i, ", ", tau.names$j, "]"), ] |> rowMeans()
    tau.mat <- matrix(tau, nrow = ncol(cov.df1), ncol = ncol(cov.df1))

    # 1. Symmetric distance component
    md <- stats::mahalanobis(cov.df1, center = mu, cov = tau.mat)
    
    # 2. Asymmetric skewness component extraction & evaluation
    alpha <- coefs[paste0("alpha[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
    X_cent <- sweep(cov.df1, 2, mu, "-")
    skew.val <- X_cent %*% alpha
    
    # 3. Combine using the Skew-Normal probability multiplier
    md.ex <- exp(beta - md/2) * 2 * stats::pnorm(skew.val)

    preds <- data.frame(cov.df[, c("x", "y")], md.ex) |> terra::rast()
    names(preds) <- paste0("Prob.", probs)
    terra::crs(preds) <- terra::crs(newdata)
  }
  
  return(preds)
}