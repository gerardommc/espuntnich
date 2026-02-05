#' @title Generate predictions from a ppmve model.
#' @description
#' Project a ppmve model onto geographic space. To do this, the user should profide
#' a fitted model, and the covariates in SpatRaster format with the same name as 
#' used to fit the model.
#' @param model = A model object of class ppmve
#' @param newdata = A SpatRaster object with covariate names contained in the covariate.names argument of ppmve function
#' @param probs = The posterior probability quantiles to be returned by predict.ppmve
#' @return Returns a single or multiple band SpatRaster object, representing point intensity as a function of distance to the estimated centroids



predict.ppmve <- function(model = NULL, newdata = NULL, probs = 0.5){
  UseMethod("predict.ppmve", model)
}

predict.ppmve.euclidean <- function(model = NULL, newdata = NULL, probs = 0.5){
  
  if(class(newdata) != "SpatRaster"){
    stop("Please provide a SpatRaster object")
  }

  if(class(newdata) == "SpatRaster"){
    cov.df <- newdata |> as.data.frame(xy = T)
    cov.df1 <- cov.df |> subset(select = c(model$call$covariates)) |> as.matrix()
  }

  if(model$call$chains == 1){
    mcl <- model$model$samples |> as.mcmc()
  }

  if(model$call$chains > 1){
    mcl <- model$model$samples |> as.mcmc.list()
    mc <- mcl[[1]]
    for(i in 2:model$call$chains){
      mc <- rbind(mc, mcl[[i]])
    }

    mcl <- mc |> as.mcmc()
  }

  if(length(probs) > 1){
    coefs <- lapply(probs, function(x){HPDinterval(mcl, prob = x)})
  }

  if(length(probs) == 1){
    coefs <- HPDinterval(mcl, prob = probs)
  }

  if(length(probs) > 1){

      p <- foreach(i = seq_along(probs), .combine = cbind) %do% {

        mu <- coefs[[i]][paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
        beta <- coefs[[i]]["beta", ] |> mean()
        tau.mat <- coefs[[i]][paste0("tau.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans() |> diag()

        md <- mahalanobis(cov.df1, center = mu, cov = tau.mat)
        md.ex <- exp(beta-md/2)
          
        return(md.ex)
          
      }

        preds <- data.frame(cov.df[, c("x", "y")], p) |> rast()
        names(preds) <- paste0("Prob.", probs)
        crs(preds) <- crs(newdata)
  }


  if(length(probs) == 1){

      mu <- coefs[paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
      beta <- coefs["beta", ] |> mean()
      tau.mat <- coefs[paste0("tau.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans() |> diag()

      md <- mahalanobis(cov.df1, center = mu, cov = tau.mat)
      md.ex <- exp(beta-md/2)

      preds <- data.frame(cov.df[, c("x", "y")], md.ex) |> rast()
      names(preds) <- paste0("Prob.", probs)
      crs(preds) <- crs(newdata)
  }

  return(preds)
}

#predict.ppmve.mahalanobis <- function(x){
#  NextMethod(x)
#}

predict.ppmve.mahalanobis <- function(model = NULL,
                                      newdata = NULL, 
                                      probs = 0.5){

  if(class(newdata) == "SpatRaster"){
    cov.df <- newdata |> as.data.frame(xy = T)
    cov.df1 <- cov.df |> subset(select = c(model$call$covariates)) |> as.matrix()
  }

  if(class(newdata) != "SpatRaster"){
    stop("Please provide a SpatRaster object")
  }


if(model$call$chains == 1){
  mcl <- model$model$samples |> as.mcmc()
}

if(model$call$chains > 1){
  mcl <- model$model$samples |> as.mcmc.list()
  mc <- mcl[[1]]
  for(i in 2:model$call$chains){
    mc <- rbind(mc, mcl[[i]])
  }

  mcl <- mc |> as.mcmc()
}

if(length(probs) > 1){
  coefs <- lapply(probs, function(x){HPDinterval(mcl, prob = x)})
}

if(length(probs) == 1){
  coefs <- HPDinterval(mcl, prob = probs)
}

tau.names <- expand.grid(i = 1:ncol(cov.df1), j = 1:ncol(cov.df1))

if(length(probs) > 1){

    p <- foreach(i = seq_along(probs), .combine = cbind) %do% {

      mu <- coefs[[i]][paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
      beta <- coefs[[i]]["beta", ] |> mean()
      tau <- coefs[[i]][paste0("tau.pres[", tau.names$i, ", ",tau.names$j, "]"), ] |> rowMeans()
      tau.mat <- matrix(tau, nrow = ncol(cov.df1), ncol = ncol(cov.df1))

      md <- mahalanobis(cov.df1, center = mu, cov = tau.mat)
      md.ex <- exp(beta-md/2)
        
      return(md.ex)
        
    }

      preds <- data.frame(cov.df[, c("x", "y")], p) |> rast()
      names(preds) <- paste0("Prob.", probs)
      crs(preds) <- crs(newdata)
  }


  if(length(probs) == 1){

    mu <- coefs[paste0("centroid.pres[", 1:ncol(cov.df1), "]"), ] |> rowMeans()
    beta <- coefs["beta", ] |> mean()
    tau <- coefs[paste0("tau.pres[", tau.names$i, ", ",tau.names$j, "]"), ] |> rowMeans()
    tau.mat <- matrix(tau, nrow = ncol(cov.df1), ncol = ncol(cov.df1))

    md <- mahalanobis(cov.df1, center = mu, cov = tau.mat)
    md.ex <- exp(beta-md/2)

    preds <- data.frame(cov.df[, c("x", "y")], md.ex) |> rast()
    names(preds) <- paste0("Prob.", probs)
    crs(preds) <- crs(newdata)
  }

  return(preds)
}
