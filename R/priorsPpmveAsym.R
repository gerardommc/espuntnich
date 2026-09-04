#' @title Generate a list of priors for espuntnich models
#' @description
#' Use a local type of CovMat ppmve or ppmveAsym model to set the prior values for a locallocal model
#' @param model A model object of class ppmve or ppmveAsym
#' @param model A model object of class ppmve or ppmveAsym
#' @param relax.precision Logical to indicate whether the precision for parameters should be reduced from what the provided model indicates
#' @param relax.factor A numeric single value or vector with length equal to the number of covariates in the model to indicate how much prior precision will be changed. 
#' Values > 1 reduce precision, and values < 1 increase it.
#' @return Returns a list with mean and precision for centroids, covariance matrices and intercept
#' @examples
#' 
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espuntnich") |> terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "pointsAsym.csv", package = "espuntnich") |> read.csv()
#' 
#' m <- ppmveAsym(points = p,
#'            covariates = r,
#'            covariate.names = names(r),
#'            CovMat = "local",
#'            Distance = "mahalanobis",
#'            no.bkgd = 5000,
#'            niter = 10000,
#'            nthin = 9,
#'            nburnin = 1000,
#'            chains = 1)
#' 
#' pri <- priors(m)
#' 
#' @export
#' @method priors ppmveAsym

priors.ppmveAsym <- function(model = NULL,
                         relax.precision = FALSE,
                         relax.factor = NULL){
  
  sts <- summary(model$model$samples)$statistics |> as.data.frame()
  
  rownms <- rownames(sts)
  
  cent.mean = sts$Mean[grep(pattern = "centroid.pres", x = rownms)]
  cent.prec = 1/(sts$SD[grep(pattern = "centroid.pres", x = rownms)])^2
  R = matrix(data = sts$Mean[grep(pattern = "tau.pres", x = rownms)],
             nrow = length(cent.mean),
             ncol = length(cent.mean),
             byrow = FALSE)
  beta.mean = sts$Mean[grep(pattern = "beta", x = rownms)]
  beta.prec = 1/(sts$SD[grep(pattern = "beta", x = rownms)])^2
  alpha.mean = sts$Mean[grep(pattern = "alpha", x = rownms)]
  alpha.prec = sts$SD[grep(pattern = "alpha", x = rownms)]

  if(relax.precision){
    if(is.null(relax.factor)){stop("Please provide a precision relaxing factor > 1")}
    cent.prec = cent.prec/relax.factor
    for(i in 1:ncol(R)){R[i, i] <- R[i, i]/relax.factor}
    beta.prec = beta.prec/relax.factor
    alpha.prec = alpha.prec/relax.factor
  }
  
  return(list(
    cent.mean = cent.mean,
    cent.prec = cent.prec,
    R = R,
    beta.mean = beta.mean,
    beta.prec = beta.prec,
    alpha.mean = alpha.mean,
    alpha.prec = alpha.prec
  ))  
}