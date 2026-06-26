#' @title Fit an asymmetric inhomogeneous Poisson Point Process Model using skew multivariate normal constraints
#' @description Fit a minimum volume ellipsoid or distance constraint with asymmetric skewness parameters as a 
#' covariate within an inhomogeneous Poisson Point Process Model to develop an Ecological Niche Model.
#' @param points A two-column data.frame with names x and y for the presence localities.
#' @param covariates A SpatRaster or imList object, containing the covariates.
#' @param covariate.names A character vector with the names of the covariates to be used in the model.
#' @param no.bkgd The number of background points to be drawn from the covariates.
#' @param bias.data The data to be used to correct observation bias.
#' @param bias.correction A character string with values "background" or "weights", indicating the bias correction method.
#' @param background.points A two-column data.frame with the x and y coordinates of background points.
#' @param samples.data A list containing three data slots: presence.data, background.data, and area.weights.
#' @param priors Informative priors list structured matching asymmetric model constraints.
#' @param CovMat A character string with values "local" or "locallocal". Relevant only for "mahalanobis" distance.
#' @param Distance A character string with values "mahalanobis" or "euclidean".
#' @param niter A numeric value indicating the number of MCMC iterations.
#' @param nburnin A numeric value indicating the number of MCMC iterations to be discarded.
#' @param nthin A numeric value indicating the interval of MCMC iterations from which posterior values will be sampled.
#' @param asCoda Logical, indicating whether the object returned by nimble run is of class coda.
#' @param chains A numeric value indicating the number of chains to be run.
#' @param WAIC Logical, whether to compute the Watanabe--Akaike information criterion.
#' @param parallel Logical, whether to run chains in parallel.
#' @param cores Numeric value indicating the number of cores to be used.
#' @param seed Numeric, indicating the random seed generator number.
#' @param remove.outer.points Logical, used to indicate if points lying outside the study window are going to be removed.
#' @param weight.bias.conf A list with configurations for the replaceQAreas function and density.ppm.
#' @return An object of class ppmveAsym and the type of distance used.
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
#' @export

ppmveAsym <- function(points = NULL,
                      covariates = NULL,
                      covariate.names = NULL,
                      no.bkgd = 5000,
                      bias.data = NULL,
                      bias.correction = NULL, 
                      background.points = NULL,
                      samples.data = NULL,
                      priors = NULL,
                      CovMat = "local",
                      Distance = "mahalanobis", 
                      niter = NULL,
                      nburnin = NULL,
                      nthin = NULL,
                      asCoda = TRUE,
                      chains = 2,
                      WAIC = TRUE,
                      parallel = FALSE,
                      cores = NULL,
                      seed = 123,
                      remove.outer.points = TRUE,
                      weight.bias.conf = list(positive = TRUE,
                                              kernel = "gaussian",
                                              sigma = NULL,
                                              varcov = NULL,
                                              weights = NULL,
                                              edge = TRUE,
                                              zo.norm = FALSE)){
  
  `%dopar%` <- foreach::`%dopar%`

  if(parallel & chains == 1){
    warning("Setting parallel as FALSE because the number of chains is 1")
    parallel <- FALSE
  }

  if(!is.null(points) & !is.null(covariates)){
    vals <- terra::extract(x = covariates, y = points, ID = FALSE)
    vals <- subset(vals, select = covariate.names)
    
    nas <- apply(vals, 2, is.na) |> apply(2, which) |> c() |> unique()
    
    if(length(nas) != 0){
      if(remove.outer.points){
        points <- points[-nas, ]
        warning(paste("Point number ", nas, " lies outside the study area and was removed from data", sep = ""))
      }

      if(!remove.outer.points){
        stop(paste("Point number ", nas, " lies outside the study area, please verify and re-run ppmveAsym", sep = ""))
      }
    }
  }

  if(!is.null(bias.correction) & !is.null(samples.data)){
    stop("Samples with data is incompatible with bias correction methods, please set either to NULL")
  }
  
  if(!is.null(background.points) & !is.null(bias.data) & !is.null(bias.correction)){
    stop("Bias correction is incompatible with user-defined background points, please set either to NULL")
  }
  
  # Preparing presence and background data
  if(!is.null(points) & !is.null(covariates)){

      names(points) <- c("x", "y")
      covariates <- covariates[[covariate.names]]
      cov.df <- as.data.frame(covariates, xy = TRUE)
      p.ext <- terra::extract(covariates, points, ID = FALSE, na.rm = TRUE) |> stats::na.omit() |> as.data.frame()
      
      iml <- espatsmo::imFromStack(covariates)
      win <- spatstat.geom::as.owin(iml[[1]])
      p.pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = win)
      
      Q <- spatstat.geom::pixelquad(p.pp)
      beg <- Q$w |> length() - iml[[1]][] |> length() + 1
      en <- Q$w |> length()

    if(is.null(bias.correction)){
      if(is.null(background.points)){
        samp <- sample(1:nrow(cov.df), no.bkgd) |> sort()
        wei <- max(Q$w)
        clim.back <- cov.df[samp, -c(1:2)]
      }

      if(!is.null(background.points)){
        if(!inherits(background.points, "matrix") & !inherits(background.points, "data.frame")){
          stop("Please provide background.points as a two-column data.frame")
        }
        if(inherits(background.points, "matrix")){
          background.points <- as.data.frame(background.points)
        }
        clim.back <- terra::extract(covariates, background.points, ID = FALSE, na.rm = TRUE) |> stats::na.omit() |> as.data.frame()
        wei <- max(Q$w)
      }
    }
      
    if(!is.null(bias.correction)){
      if (bias.correction == "weights") {
        if (inherits(bias.data, "data.frame")) {
          
          bias.ppp <- spatstat.geom::ppp(x = bias.data$x, 
                                         y = bias.data$y, window = win)
          
          dens.r <- terra::rast(spatstat.geom::density.ppp(bias.ppp, 
                                                           positive = weight.bias.conf$positive, kernel = weight.bias.conf$kernel, 
                                                           sigma = weight.bias.conf$sigma, varcov = weight.bias.conf$varcov, 
                                                           weights = weight.bias.conf$weights, edge = weight.bias.conf$edge))
          
          bias.df <- terra::as.data.frame(dens.r, xy = TRUE)
          ux <- bias.df$x |> unique() |> sort()
          uy <- bias.df$y |> unique() |> sort()
          
          dx <- ux[2] - ux[1]
          dy <- uy[2] - uy[1]
          
          area <- dx * dy
          
          sum.bias <- terra::global(bias.data, sum, na.rm = TRUE)
          
          ones <- terra::classify(bias.data, rcl = matrix(c(-Inf, Inf, 1), ncol = 3))
          
          sum.ones <- terra::global(ones, sum, na.rm = TRUE)
          
          bias.data <- bias.data / sum.bias$sum
          p.area <- ones / sum.ones$sum
          
          Ratio <- bias.data / p.area
          
          area.weights <- Ratio * area
          
          samp <- sort(sample(1:nrow(cov.df), no.bkgd))
          wei <- terra::extract(area.weights, cov.df[samp, c("x", "y")], ID = FALSE)
          
          nas <- is.na(wei[, 1])
          samp <- samp[!nas]
          wei <- wei[!nas, 1]
        }
        if (inherits(bias.data, "SpatRaster")) {

          bias.df <- terra::as.data.frame(bias.data, xy = TRUE)
          ux <- bias.df$x |> unique() |> sort()
          uy <- bias.df$y |> unique() |> sort()
          
          dx <- ux[2] - ux[1]
          dy <- uy[2] - uy[1]
          
          area <- dx * dy
          
          sum.bias <- terra::global(bias.data, sum, na.rm = TRUE)
          
          ones <- terra::classify(bias.data, rcl = matrix(c(-Inf, Inf, 1), ncol = 3))
          
          sum.ones <- terra::global(ones, sum, na.rm = TRUE)
          
          bias.data <- bias.data / sum.bias$sum
          p.area <- ones / sum.ones$sum
          
          Ratio <- bias.data / p.area
          
          area.weights <- Ratio * area

          samp <- sort(sample(1:nrow(cov.df), no.bkgd))
          wei <- terra::extract(area.weights, cov.df[samp, c("x", "y")], ID = FALSE)
          
          nas <- is.na(wei[, 1])
          samp <- samp[!nas]
          wei <- wei[!nas, 1]
        }
      }
        
      if(bias.correction == "background"){
        if(inherits(bias.data, "data.frame")){
            bias.ppp <- spatstat.geom::ppp(x = bias.data$x, y = bias.data$y, window = win)
            dens.r <- spatstat.geom::density.ppp(bias.ppp,
                                                 positive = weight.bias.conf$positive,
                                                 kernel = weight.bias.conf$kernel,
                                                 sigma = weight.bias.conf$sigma,
                                                 varcov = weight.bias.conf$varcov,
                                                 weights = weight.bias.conf$weights,
                                                 edge = weight.bias.conf$edge) |> terra::rast()
            
            dens.df <- as.data.frame(dens.r, xy = TRUE)
            area.weights <- iml[[1]]
            area.weights[] <- Q$w[beg:en]
            weights.r <- area.weights |> terra::rast()
            samp <- sample(1:nrow(dens.df), no.bkgd, prob = dens.df[, 3]) |> sort()
            wei <- max(Q$w)
        }
          
        if(inherits(bias.data, "SpatRaster")){
            bias.data <- terra::resample(bias.data, covariates[[1]]) |> espatsmo::ZeroOneNorm()
            bias.df <- as.data.frame(bias.data, xy = TRUE)
            area.weights <- iml[[1]]
            area.weights[] <- Q$w[beg:en]
            weights.r <- area.weights |> terra::rast()
            samp <- sample(1:nrow(bias.df), no.bkgd, prob = bias.df[, 3]) |> sort()
            wei <- max(Q$w)
        }
      }
      clim.back <- cov.df[samp, -c(1:2)]
    }
  }
    
  if(!is.null(samples.data$presence.data) & !is.null(samples.data$background.data)){
      clim.back <- samples.data$background.data[, covariate.names]
      wei <- samples.data$area.weights
      p.ext <- samples.data$presence.data[, covariate.names]
  }

  # Configuring data, constants and parameters for asymmetric versions
  if(Distance == "mahalanobis"){
    
    # Asymmetric Mahalanobis monitors track skewness parameters (alpha)
    parms <- c("centroid.pres", "mu.back", "tau.pres", "alpha", "beta")
    
    if(CovMat == "local"){
      if(is.null(priors)){
        constants <- list(n.clim = ncol(clim.back),
                          R = diag(ncol(clim.back)),
                          n.data = nrow(points) + nrow(clim.back),
                          n.pres = nrow(points),
                          cent.mean = rep(0, ncol(clim.back)),
                          cent.prec = rep(0.1, ncol(clim.back)),
                          mu.b.mean = rep(0, ncol(clim.back)),
                          mu.b.prec = rep(1.0E-4, ncol(clim.back)),
                          alpha.mean = rep(0, ncol(clim.back)),
                          alpha.prec = rep(0.1, ncol(clim.back)),
                          beta.mean = 0,
                          beta.prec = 1.0E-4)
      } else {
        constants <- list(n.clim = ncol(clim.back),
                          R = priors$R,
                          n.data = nrow(points) + nrow(clim.back),
                          n.pres = nrow(points),
                          cent.mean = priors$cent.mean,
                          cent.prec = priors$cent.prec,
                          mu.b.mean = priors$mu.b.mean,
                          mu.b.prec = priors$mu.b.prec,
                          alpha.mean = priors$alpha.mean,
                          alpha.prec = priors$alpha.prec,
                          beta.mean = priors$beta.mean,
                          beta.prec = priors$beta.prec)            
      }
    }
    
    if(CovMat == "locallocal"){
      parms <- c("centroid.pres", "tau.pres", "alpha", "beta")
      if(is.null(priors)){
        constants <- list(n.clim = ncol(clim.back),
                          R = diag(ncol(clim.back)),
                          n.data = nrow(points) + nrow(clim.back),
                          cent.mean = rep(0, ncol(clim.back)),
                          cent.prec = rep(0.1, ncol(clim.back)),
                          alpha.mean = rep(0, ncol(clim.back)),
                          alpha.prec = rep(0.1, ncol(clim.back)),
                          beta.mean = 0,
                          beta.prec = 1.0E-4)
      } else {
        constants <- list(n.clim = ncol(clim.back),
                          R = priors$R,
                          n.data = nrow(points) + nrow(clim.back),
                          cent.mean = priors$cent.mean,
                          cent.prec = priors$cent.prec,
                          alpha.mean = priors$alpha.mean,
                          alpha.prec = priors$alpha.prec,
                          beta.mean = priors$beta.mean,
                          beta.prec = priors$beta.prec)
      }
    }

    if(length(wei) == 1){
      constants$w <- c(rep(1/wei, nrow(points)), rep(wei, nrow(clim.back)))
    } else {
      constants$w <- c(rep(1/(stats::median(wei)), nrow(points)), wei)
    }
              
    inits <- list(centroid.pres = constants$cent.mean,
                  tau.pres = constants$R,
                  alpha = constants$alpha.mean,
                  beta = constants$beta.mean)
    if(CovMat != "locallocal") inits$mu.back <- constants$mu.b.mean
      
    data <- list(lambda = c(rep(1L, nrow(points)), rep(0L, nrow(clim.back))),
                 clim = rbind(p.ext, clim.back))
  }
  
  if(Distance == "euclidean"){
    
    # Asymmetric Euclidean tracks shape metrics (alpha) alongside base metrics
    parms <- c("centroid.pres", "tau.pres", "alpha", "beta")
    
    if(is.null(priors)){
        constants <- list(n.clim = ncol(clim.back),
                          n.data = nrow(points) + nrow(clim.back),
                          cent.mean = rep(0, ncol(clim.back)),
                          cent.prec = rep(0.1, ncol(clim.back)),
                          tau.min = rep(0, ncol(clim.back)),
                          tau.max = rep(100, ncol(clim.back)),
                          alpha.mean = rep(0, ncol(clim.back)),
                          alpha.prec = rep(0.1, ncol(clim.back)),
                          beta.mean = 0,
                          beta.prec = 1.0E-4)
    } else {
        constants <- list(n.clim = ncol(clim.back),
                          n.data = nrow(points) + nrow(clim.back),
                          cent.mean = priors$cent.mean,
                          cent.prec = priors$cent.prec,
                          tau.min = priors$tau.min,
                          tau.max = priors$tau.max,
                          alpha.mean = priors$alpha.mean,
                          alpha.prec = priors$alpha.prec,
                          beta.mean = priors$beta.mean,
                          beta.prec = priors$beta.prec)
    } 
    
    if(length(wei) == 1){
      constants$w <- c(rep(1/wei, nrow(points)), rep(wei, nrow(clim.back)))
    } else {
      constants$w <- c(rep(1/(stats::median(wei)), nrow(points)), wei)
    }
    
    inits <- list(centroid.pres = constants$cent.mean,
                  tau.pres = rep(1, constants$n.clim),
                  alpha = constants$alpha.mean,
                  beta = constants$beta.mean)
    
    data <- list(lambda = c(rep(1L, nrow(points)), rep(0L, nrow(clim.back))),
                 clim = rbind(p.ext, clim.back))
  }
  
  # Selecting the specified Nimble Asymmetric Model
  if(Distance == "mahalanobis"){
    if(CovMat == "local"){
      modelCode <- LocalMahalAsym
    }
    if(CovMat == "locallocal"){
      modelCode <- LocalLocalMahalAsym
    }
  }
  
  if(Distance == "euclidean"){
    modelCode <- EuclidAsym
  }
  
  # Coordinate dimensions compilation setup
  dim_list <- if(Distance == "mahalanobis") {
    list(covMat.pres = c(constants$n.clim, constants$n.clim))
  } else {
    list()
  }

  model <- nimble::nimbleModel(code = modelCode,
                               constants = constants,
                               data = data,
                               inits = inits,
                               dimensions = dim_list)
  
  cmodel <- nimble::compileNimble(model)
  conf.model <- nimble::configureMCMC(model, monitors = parms, enableWAIC = TRUE)
  MCMC <- nimble::buildMCMC(conf.model)
  cMCMC <- nimble::compileNimble(MCMC, project = cmodel)
  
  set.seed(seed)
  
  if(parallel){
    doParallel::registerDoParallel(cores = cores)
    
    run.list <- foreach::foreach(i = 1:chains, .packages = "nimble") %dopar% {
      ch <- nimble::runMCMC(cMCMC, niter = niter, nburnin = nburnin, thin = nthin, samplesAsCodaMCMC = asCoda, nchains = 1, WAIC = WAIC)
      return(ch)
    }
    
    mc.list <- list()
    for(i in seq_along(run.list)){
      temp <- run.list[[i]]
      mc.list[[i]] <- temp$samples
    }
    names(mc.list) <- paste0("chain", seq_along(run.list))
    
    run <- run.list[[1]]
    run$samples <- coda::as.mcmc.list(mc.list)

    if(!is.null(background.points)){
      num.back <- nrow(background.points)
      b.points <- background.points
    } else if(!is.null(samples.data)){
      num.back <- nrow(samples.data$background.data)
      b.points <- NULL
    } else {
      num.back <- no.bkgd
      b.points <- cov.df[samp, c(1:2)]
    }

    ret.list <- list(model = run,
                     call = list(Distance = Distance,
                                 CovMat = CovMat,
                                 no.bkgd = num.back,
                                 covariates = covariate.names,
                                 niter = niter,
                                 nburnin = nburnin,
                                 nthin = nthin,
                                 asCoda = asCoda,
                                 chains = chains,
                                 WAIC = WAIC,
                                 parallel = parallel,
                                 cores = cores,
                                 seed = seed),
                     bkgd.points = b.points)
    
    if(Distance == "mahalanobis"){
      class(ret.list) <- c("ppmve", Distance, CovMat)
    }
    
    if(Distance == "euclidean"){
      class(ret.list) <- c("ppmve", Distance)
    }
    return(ret.list)
    
  } else {
    run <- nimble::runMCMC(cMCMC, niter = niter, nburnin = nburnin, thin = nthin, samplesAsCodaMCMC = asCoda, nchains = chains, WAIC = WAIC)

    if(!is.null(background.points)){
      num.back <- nrow(background.points)
      b.points <- background.points
    } else if(!is.null(samples.data)){
      num.back <- nrow(samples.data$background.data)
      b.points <- NULL
    } else {
      num.back <- no.bkgd
      b.points <- cov.df[samp, c(1:2)]
    }

    ret.list <- list(model = run,
                     call = list(Distance = Distance,
                                 CovMat = CovMat,
                                 no.bkgd = num.back,
                                 covariates = covariate.names,
                                 niter = niter,
                                 nburnin = nburnin,
                                 nthin = nthin,
                                 asCoda = asCoda,
                                 chains = chains,
                                 WAIC = WAIC,
                                 parallel = parallel,
                                 cores = cores,
                                 seed = seed),
                     bkgd.points = b.points)
    
    if(Distance == "mahalanobis"){
      class(ret.list) <- c("ppmve", Distance, CovMat)
    }
    
    if(Distance == "euclidean"){
      class(ret.list) <- c("ppmve", Distance)
    }
    return(ret.list)
  }
}