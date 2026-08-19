#' @title Fit an inhomogeneous Poisson Point Process Model as a function of distance to an ellipsoid's centroid
#' @description Fit a minimum volume elipsoid as a covariate within an inhomogeneous Poisson Point Process Model 
#' to develop an Ecological Niche Model using point locations and at least two raster covariates.
#' @param points A two-column data.frame with names x and y for the presence localities.
#' @param covariates A SpatRaster or imList object, containing the covariates.
#' @param covariate.names A character vector with the names of the covariaes to be used in the model.
#' @param no.bkgd The number of background points to be drawn from the covariates.
#' @param bias.data The data to be used to correct observation bias. This data could be a data frame of fieldwork sampling localities, or a raster layer describing the variability of observation effort.
#' @param bias.correction A character string with values "background" or "weights", indicating the bias correction method. For "background", bias is done by selecting sampling localities from more intensely sampled areas. For weights, the size of the areas are modified according to the ratio of uniform vs. realised sampling effort.
#' @param background.points A two-column data.frame with the x and y coordinates of background or quadrature points used to sample covariate space.
#' @param samples.data A list containing three data slots: (1) presence.data, a data.frame where each column contains the environmental conditions of each presence point, and where column names
#' must contain the names supplied in covariate names; (2) background.data,a data.frame where each column contains the environmental conditions of background locations, and where 
#' column names must contain the supplied covariate names; (3) area.weights, a single numeric value or vector containing the size of covariate pixels. For example: list(presence.data = NULL,
#' background.data = NULL, area.weights = NULL)
#' @param priors By default it is NULL, and flat priors are assumed for all centroid coordinates or precision components. If more informative priors are desired, a list with data slots
#' named after each parameter estimated by the supported models must be provided. The variations for mahalanobis distance models are (with parameters which must be contained in the list 
#' between brackets): (1) local (centroid.pres, mu.back, tau.pres, beta); (2) locallocal (centroid.pres, tau.pres, beta); 
#' For the euclidean model, the estimated parameters are "centroid.pres", "tau.pres" and "beta". "centrod.pres" is a real-valued vector equal to those used in teh mahalanobis models; "tau.pres" is a 
#' positive real-valued vector of length equal to the length of "covariate.names", and "beta" is the same as in the mahalanobis models.
#' @param CovMat A character string with values "local", "locallocal", to configure whether the covariance matrix is parameterised from the precence localities or point process only. This argument is only relevant for "mahalanobis" distance.
#' @param Distance A character string with values"mahalanobis" or "euclidean", indicating the type of distance to be calculated.
#' @param niter A numeric value indicating the number of MCMC iterations.
#' @param nburnin A numeric value indicating the number of MCMC iterations to be discarded at the beginning of each chain.
#' @param nthin A numeric value indicating the interval of MCMC iterations from which posterior values will be sampled
#' @param asCoda Logical, indicating whether the objecct returned by nimble rum is of class coda. This shoule be set to T to use predict methods.
#' @param chains A numeric value indicating the number of chains to be run. If convergence diagnostics are to be run, the number should be at least 2.
#' @param WAIC Logical, whether to compute the Watanabe--Akaike information criterion, T.
#' @param parallel Logical, whether to run chains in parallel, F.
#' @param cores Numeric value indicating the number of cores to be used for running the chains. Should be equal to "chains".
#' @param seed Numeric, idincating the random seed generator number, 123.
#' @param remove.outer.points Logical, used to indicate iff points ling outside the study window are going to be removed from model run.
#' @param weight.bias.conf = A list with default entries, positive = TRUE, kernel = "gaussian", sigma = NULL, varcov = NULL, weights = NULL, edge = TRUE, which are used to configure the replaceQAreas function and density.ppm, only relevant if bias.correction = "weights" and the class of bias.data is data.frame
#' @return An object of class ppmve and the type of distance and covariance matrix used.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
#' 
#' # Model with random sampling
#' m <- ppmve(points = pr,
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
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#' 
#' sb <- system.file("extdata", "BiasSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pb <- p[sb$Samples, ]
#' 
#' # Model with biased sampling
#' mb <- ppmve(points = pb,
#'            covariates = r,
#'            covariate.names = names(r),
#'            CovMat = "local",
#'            Distance = "mahalanobis",
#'            bias.data = bias,
#'            bias.correction = "weights",
#'            no.bkgd = 5000,
#'            niter = 10000,
#'            nthin = 9,
#'            nburnin = 1000,
#'            chains = 1)
#' @export

ppmve <- function(points = NULL,
                  covariates = NULL,
                  covariate.names = NULL,
                  no.bkgd = 5000,
                  bias.data = NULL,
                  bias.correction = NULL, #options = "background", "weights"
                  background.points = NULL,
                  samples.data = NULL,
                  priors = NULL,
                  CovMat = "local",
                  Distance = "mahalanobis", #options = "euclidean"
                  niter = NULL,
                  nburnin = NULL,
                  nthin = NULL,
                  asCoda = TRUE,
                  chains = 2,
                  WAIC = T,
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
        stop(paste("Point number ", nas, " lies outside the study area, please verify and re-run ppmve", sep = ""))
      }
    }
  }


  if(!is.null(bias.correction) & !is.null(samples.data)){
    stop("Samples with data is incompatible with bias correction methods, please set either to NULL")
  }
  
  if(!is.null(background.points) & !is.null(bias.data) & !is.null(bias.correction)){
    stop("Bias correction is incompatible with user-defined background points, please set either to NULL")
  }
  
#Preparing presence and background data
  if(!is.null(points) & !is.null(covariates)){

      names(points) <- c("x", "y")
      
      covariates <- covariates[[covariate.names]]
      
      cov.df <- as.data.frame(covariates, xy = TRUE)
      p.ext <- terra::extract(covariates, points, ID = F, na.rm = TRUE) |> stats::na.omit() |> as.data.frame()
      
      iml <- espatsmo::imFromStack(covariates)
      win <- spatstat.geom::as.owin(iml[[1]])
      p.pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = win)
      
      #Calculating weights
      Q <- spatstat.geom::pixelquad(p.pp)
      
      beg <- Q$w |> length() - iml[[1]][] |>length() + 1
      en <- Q$w |> length()

    if(is.null(bias.correction)){
      if(is.null(background.points)){
        samp <- sample(1:nrow(cov.df), no.bkgd) |> sort()
        wei <- max(Q$w) * nrow(cov.df)/length(samp)
        clim.back <- cov.df[samp, -c(1:2)]
      }

      if(!is.null(background.points)){

        if(!inherits(background.points, "matrix") & !inherits(background.points, "data.frame")){
          stop("Please provide background.points as a two-column data.frame")
        }

        if(inherits(background.points, "matrix")){
          background.points <- as.data.frame(background.points)
        }

        clim.back <- terra::extract(covariates, background.points, ID = F, na.rm = T) |> stats::na.omit() |> as.data.frame()
        wei <- max(Q$w) * nrow(cov.df)/nrow(background.points)
      }
    }
      
    if(!is.null(bias.correction)){
        
        #Background data with bias correction
        #Bias correction based on are weights
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
          wei <- wei[!nas, 1] * nrow(cov.df)/length(samp)
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
          wei <- wei[!nas, 1] * nrow(cov.df)/length(samp)
        }
      }
        
        #Bias correction based on location of  background data
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
            
            wei <- max(Q$w) * nrow(cov.df)/length(samp)
        }
          
        if(inherits(bias.data, "SpatRaster")){
            bias.data <- terra::resample(bias.data, covariates[[1]]) |> espatsmo::ZeroOneNorm()
            
            bias.df <- as.data.frame(bias.data, xy = TRUE)
            
            area.weights <- iml[[1]]
            area.weights[] <- Q$w[beg:en]
            weights.r <- area.weights |> terra::rast()
            
            samp <- sample(1:nrow(bias.df), no.bkgd, prob = bias.df[, 3]) |> sort()
            
            wei <- max(Q$w) * nrow(cov.df)/length(samp)
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

  # Configuring data, constants and parameters for each model type 
  
  if(Distance == "mahalanobis"){

    if(CovMat == "local"){

      parms <- c("centroid.pres",
                "mu.back",
                "tau.pres",
                "beta")
      
      if(is.null(priors)){
          constants <- list(n.clim = ncol(clim.back),
                R = diag(ncol(clim.back)),
                n.data = nrow(points) + nrow(clim.back),
                n.pres = nrow(points),
                cent.mean = rep(0, ncol(clim.back)),
                cent.prec = rep(0.1, ncol(clim.back)),
                mu.b.mean = rep(0, ncol(clim.back)),
                mu.b.prec = rep(1.0E-4, ncol(clim.back)),
                beta.mean = 0,
                beta.prec = 1.0E-4)
        }

        if(!is.null(priors)){
          constants <- list(n.clim = ncol(clim.back),
                R = priors$R,
                n.data = nrow(points) + nrow(clim.back),
                n.pres = nrow(points),
                cent.mean = priors$cent.mean,
                cent.prec = priors$cent.prec,
                mu.b.mean = priors$mu.b.mean,
                mu.b.prec = priors$mu.b.prec,
                beta.mean = priors$beta.mean,
                beta.prec = priors$beta.prec)            
        }

      if(length(wei) == 1){
          constants$w <- c(rep(1/wei, nrow(points)), rep(wei, nrow(clim.back)))
      }
      
      if(length(wei) > 1){
          constants$w <- c(rep(1/(stats::median(wei)), nrow(points)), wei)
      }
      
      inits <- list(centroid.pres = constants$cent.mean,
                    mu.back = constants$mu.b.mean,
                    tau.pres = constants$R,
                    beta = constants$beta.mean)
      
      data <- list(lambda = c(rep(1L, nrow(points)), rep(0L, nrow(clim.back))),
                  clim = rbind(p.ext, clim.back))      
    }
    
    if(CovMat == "locallocal"){
      
      parms <- c("centroid.pres",
                "tau.pres",
                "beta")
      
      if(is.null(priors)){
        constants <- list(n.clim = ncol(clim.back),
              R = diag(ncol(clim.back)),
              n.data = nrow(points) + nrow(clim.back),
              cent.mean = rep(0, ncol(clim.back)),
              cent.prec = rep(0.1, ncol(clim.back)),
              beta.mean = 0,
              beta.prec = 1.0E-4)
      }

      if(!is.null(priors)){
        constants <- list(n.clim = ncol(clim.back),
              R = priors$R,
              n.data = nrow(points) + nrow(clim.back),
              cent.mean = priors$cent.mean,
              cent.prec = priors$cent.prec,
              beta.mean = priors$beta.mean,
              beta.prec = priors$beta.prec)
      }

      if(length(wei) == 1){
          constants$w <- c(rep(1/wei, nrow(points)), rep(wei, nrow(clim.back)))
      }
      
      if(length(wei) > 1){
          constants$w <- c(rep(1/(stats::median(wei)), nrow(points)), wei)
      }
      
      inits <- list(centroid.pres = constants$cent.mean,
                    tau.pres = constants$R,
                    beta = constants$beta.mean)
      
      data <- list(lambda = c(rep(1L, nrow(points)), rep(0L, nrow(clim.back))),
                  clim = rbind(p.ext, clim.back))
      
    }
  }
  
  if(Distance == "euclidean"){
    
    parms <- c("centroid.pres",
              "tau.pres",
              "beta")
    
    if(is.null(priors)){
        constants <- list(n.clim = ncol(clim.back),
                          n.data = nrow(points) + nrow(clim.back),
                          cent.mean = rep(0, ncol(clim.back)),
                          cent.prec = rep(0.1, ncol(clim.back)),
                          tau.min = rep(0, ncol(clim.back)),
                          tau.max = rep(100, ncol(clim.back)),
                          beta.mean = 0,
                          beta.prec = 1.0E-4)
      }

      if(!is.null(priors)){
        constants <- list(n.clim = ncol(clim.back),
                          n.data = nrow(points) + nrow(clim.back),
                          cent.mean = priors$cent.mean,
                          cent.prec = priors$cent.prec,
                          tau.min = priors$tau.min,
                          tau.max = priors$tau.max,
                          beta.mean = priors$beta.mean,
                          beta.prec = priors$beta.prec)
      } 
    
    if(length(wei) == 1){
      constants$w <- c(rep(1/wei, nrow(points)), rep(wei, nrow(clim.back)))
    }
    
    if(length(wei) > 1){
      constants$w <- c(rep(1/(stats::median(wei)), nrow(points)), wei)
    }
    
    inits <- list(centroid.pres = constants$cent.mean,
                  tau.pres = rep(1, constants$n.clim),
                  beta = constants$beta.mean)
    
    data <- list(lambda = c(rep(1L, nrow(points)), rep(0L, nrow(clim.back))),
                clim = rbind(p.ext, clim.back))
  }
  
  #Loading the specified Nimble model
  if(Distance == "mahalanobis"){
    
    if(CovMat == "local"){
      modelCode <- LocalMahal
    }
    
    if(CovMat == "locallocal"){
      modelCode <- LocalLocalMahal
    }
  }
  
  if(Distance == "euclidean"){
    modelCode <- Euclid
  }
  
  model <- nimble::nimbleModel(code = modelCode,
                               constants = constants,
                               data = data,
                               inits = inits,
                               dimensions = list(covMat.pres = c(constants$n.clim, constants$n.clim)))
  
  cmodel <- nimble::compileNimble(model)
  
  conf.model <- nimble::configureMCMC(model, monitors = parms, enableWAIC = T)
  
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
    }

    if(!is.null(samples.data)){
      num.back <- nrow(samples.data$background.data)
      b.points <- NULL
    }

    if(is.null(background.points) & is.null(samples.data)){
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
    }

    if(!is.null(samples.data)){
      num.back <- nrow(samples.data$background.data)
      b.points <- NULL
    }

    if(is.null(background.points) & is.null(samples.data)){
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

