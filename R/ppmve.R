#' @title Fit an inhomogeneous Poisson Point Process Model as a function of distance to an ellipsoid's centroid
#' @description
#' @param points = A two-column data.frame with names x and y for the presence localities.
#' @param covariates = A SpatRaster or imList object, containing the covariates.
#' @param covariate.names = A character vector with the names of the covariaes to be used in the model.
#' @param no.bkgd = The number of background points to be drawn from the covariates.
#' @param bias.data = The data to be used to correct observation bias. This data could be a data frame of fieldwork sampling localities, or a raster layer describing the variability of observation effort.
#' @param bias.correction = A character string with values "background" or "weights", indicating the bias correction method. For "background", bias is done by selecting sampling localities from more intensely sampled areas. For weights, the size of the areas are modified according to the ratio of uniform vs. realised sampling effort.
#' @param CovMat = A character string with values "local" or "global", to configure whether the covariance matrix is parameterised from the point process only or separately from values at sampling localities. This argument is only relevant for "mahalanobis" distance.
#' @param Distance =  A character string with values"mahalanobis" or "euclidean", indicating the type of distance to be calculated.
#' @param niter = A numeric value indicating the number of MCMC iterations.
#' @param nburnin = A numeric value indicating the number of MCMC iterations to be discarded at the beginning of each chain.
#' @param thin = A numeric value indicating the interval of MCMC iterations from which posterior values will be sampled
#' @param asCoda = Logical, indicating whether the objecct returned by nimble rum is of class coda. This shoule be set to T to use predict methods.
#' @param chains = A numeric value indicating the number of chains to be run. If convergence diagnostics are to be run, the number should be at least 2.
#' @param WAIC = Logical, whether to compute the Watanabe--Akaike information criterion, T.
#' @param parallel = Logical, whether to run chains in parallel, F.
#' @param cores = Numeric value indicating the number of cores to be used for running the chains. Should be equal to "chains".
#' @param seed = Numeric, idincating the random seed generator number, 123.
#' @param weight.bias.conf = A list with default entries, nsim =  39, positive = TRUE, kernel = "gaussian", sigma = NULL, varcov = NULL, weights = NULL, edge = TRUE, which are used to configure the replaceQAreas function and density.ppm, only relevant if bias.correction = "weights" and the class of bias.data is data.frame
#' @return An object of class ppmve and the type of distance and covariance matrix used.

ppmve <- function(points = NULL,
                  covariates = NULL,
                  covariate.names = NULL,
                  no.bkgd = 5000,
                  bias.data = NULL,
                  bias.correction = NULL, #options = "background", "weights"
                  CovMat = "local",
                  Distance = "mahalanobis", #options = "euclidean"
                  niter = NULL,
                  nburnin = NULL,
                  thin = NULL,
                  asCoda = T,
                  chains = 2,
                  WAIC = T,
                  parallel = F,
                  cores = NULL,
                  seed = 123,
                  weight.bias.conf = list(
                    nsim =  39,
                    positive = TRUE,
                    kernel = "gaussian",
                    sigma = NULL,
                    varcov = NULL,
                    weights = NULL,
                    edge = TRUE)){


  if(is.null(points) | is.null(covariates) | is.null(covariate.names)){
    stop("Please specify valid inputs for points, covariates and covariate.names")
  }

  names(points) <- c("x", "y")

  covariates <- covariates[[covariate.names]]

  cov.df <- as.data.frame(covariates, xy = T)
  p.ext <- terra::extract(covariates, points, ID = F, na.rm = T) |> na.omit() |> as.data.frame()

  iml <- espatsmo::imFromStack(covariates)
  win <- spatstat.geom::as.owin(iml[[1]])
  p.pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = win)

#Calculating weights
  Q <- spatstat.geom::pixelquad(p.pp)

  beg <- Q$w |> length() - iml[[1]][] |>length() + 1
  en <- Q$w |> length()

  if(is.null(bias.correction)){
    samp <- sample(1:nrow(cov.df), no.bkgd) |> sort()
    wei <- max(Q$w)
  }

  if(!is.null(bias.correction)){

    #Background data with bias correction
    #Bias correction based on are weights
    if(bias.correction == "weights"){
      if(class(bias.data) == "data.frame"){
        Qa <- espatsmo::replaceQAreas(Q = Q,
                            bias.data = bias.data,
                            im = iml[[1]],
                            nsim =  weight.bias.conf$nsim,
                            positive = weight.bias.conf$positive,
                            kernel = weight.bias.conf$kernel,
                            sigma = weight.bias.conf$sigma,
                            varcov = weight.bias.conf$sigma,
                            weights = weight.bias.conf$weights,
                            edge = weight.bias.conf$edge)

        area.weights <- iml[[1]]
        area.weights[] <- Qa$w[beg:en]
        weights.r <- area.weights |> terra::rast()

        samp <- sample(1:nrow(cov.df), no.bkgd) |> sort()

        wei <- terra::extract(weights.r, cov.df[samp, c("x", "y")])[,2]

        nas <- wei |> is.na()

        samp <- samp[!nas]

        wei <- wei[!nas]
    }

  if(class(bias.data) == "SpatRaster"){

    bias.data <- terra::resample(bias.data, covariates[[1]]) |> espatsmo::ZeroOneNorm()
    bias.df <- as.data.frame(bias.data, xy = TRUE)
    ids.bias <- sample(1:nrow(bias.df), no.bkgd, prob = bias.df[, 3]) |> sort()
    locs.bias <- bias.df[ids.bias, ]

    Qa <- espatsmo::replaceQAreas(Q = Q,
                        bias.data = locs.bias,
                        im = iml[[1]],
                        nsim =  weight.bias.conf$nsim,
                        positive = weight.bias.conf$positive,
                        kernel = weight.bias.conf$kernel,
                        sigma = weight.bias.conf$sigma,
                        varcov = weight.bias.conf$sigma,
                        weights = weight.bias.conf$weights,
                        edge = weight.bias.conf$edge)

    area.weights <- iml[[1]]
    area.weights[] <- Qa$w[beg:en]
    weights.r <- area.weights |> terra::rast()
    samp <- sample(1:nrow(cov.df), no.bkgd) |> sort()
    wei <- terra::extract(weights.r, cov.df[samp, c("x", "y")])
    nas <- wei |> is.na()
    samp <- samp[!nas]
    wei <- wei[!nas]
    }
  }

#Bias correction based on location of  background data
  if(bias.correction == "background"){
    if(class(bias.data) == "data.frame"){
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

    if(class(bias.data) == "SpatRaster"){
      bias.data <- terra::resample(bias.data, covariates[[1]]) |> espatsmo::ZeroOneNorm()

      bias.df <- as.data.frame(bias.data, xy = TRUE)

      area.weights <- iml[[1]]
      area.weights[] <- Q$w[beg:en]
      weights.r <- area.weights |> terra::rast()

      samp <- sample(1:nrow(bias.df), no.bkgd, prob = bias.df[, 3]) |> sort()

      wei <- max(Q$w)
      }
    }
  }

  clim.back <- cov.df[samp, -c(1:2)]

  # Loading contstants for each model type
  if(Distance == "mahalanobis"){
    if(length(wei) == 1){
      constants <- list(n.clim = ncol(clim.back),
                        w = c(rep(1/wei, nrow(points)),
                              rep(wei, nrow(clim.back))),
                        R = diag(ncol(clim.back)),
                        n.data = nrow(points) + nrow(clim.back),
                        n.pres = nrow(points),
                        n.back = nrow(clim.back))
    } else {
      constants <- list(n.clim = ncol(clim.back),
                        w = c(rep(1/(median(wei)), nrow(points)), wei),
                        R = diag(ncol(clim.back)),
                        n.data = nrow(points) + nrow(clim.back),
                        n.pres = nrow(points),
                        n.back = nrow(clim.back))
    }
  }

  if(Distance == "euclidean"){
    if(length(wei) == 1){
      constants <- list(n.clim = ncol(clim.back),
                        w = c(rep(1/wei, nrow(points)),
                              rep(wei, nrow(clim.back))),
                        n.data = nrow(points) + nrow(clim.back),
                        n.back = nrow(clim.back))
    } else {
      constants <- list(n.clim = ncol(clim.back),
                        w = c(rep(1/(median(wei)), nrow(points)), wei),
                        n.data = nrow(points) + nrow(clim.back),
                        n.back = nrow(clim.back))
    }
  }

  # Configure parms, data and inits for model types

  if(Distance == "mahalanobis"){
      parms <- c("centroid.pres",
                  "mu.back",
                  "tau.pres",
                  "beta")

      inits <- list(centroid.pres = rep(0, constants$n.clim),
              mu.back = rep(0, constants$n.clim),
              tau.pres = diag(1, nrow = constants$n.clim),
              beta = 0)
    }

  if(Distance == "euclidean"){
    inits <- list(centroid.pres = rep(0, constants$n.clim),
              tau.pres = rep(1, constants$n.clim),
              beta = 0)

    parms <- c("centroid.pres",
           "tau.pres",
           "beta")
  }


  data <- list(lambda = c(rep(1L, nrow(points)), rep(0L, nrow(clim.back))),
               clim = rbind(p.ext, clim.back),
               n.clim = ncol(cov.df[, -(1:2)]))

  #Loading the specified Nimble model
  if(Distance == "mahalanobis"){
    if(CovMat == "global"){
      modelCode <- GlobalMahal
    }

    if(CovMat == "local"){
      modelCode <- GlobalMahal
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
    run <- doParallel::foreach(i = 1:chains) %dopar% {
      nimble::runMCMC(cMCMC, niter = niter, nburnin = nburnin, thin = thin, samplesAsCodaMCMC = asCoda, nchains = 1, WAIC = F)
    }

    ret.list <- list(model = run,
                     call = list(Distance = Distance,
                                CovMat = CovMat,
                                no.bkgd = no.bkgd,
                                covariates = names(covariates),
                                niter = niter,
                                nburnin = nburnin,
                                thin = thin,
                                asCoda = asCoda,
                                chains = chains,
                                WAIC = WAIC,
                                parallel = parallel,
                                cores = cores,
                                seed = seed),
                                bkgd.points = cov.df[samp, c(1:2)])

    class(ret.list) <- c("ppmve", Distance, CovMat)

    return(ret.list)
  } else {
    run <- nimble::runMCMC(cMCMC, niter = niter, nburnin = nburnin, thin = thin, samplesAsCodaMCMC = asCoda, nchains = chains, WAIC = WAIC)
    ret.list <- list(model = run,
                     call = list(Distance = Distance,
                                CovMat = CovMat,
                                no.bkgd = no.bkgd,
                                covariates = covariate.names,
                                niter = niter,
                                nburnin = nburnin,
                                thin = thin,
                                asCoda = asCoda,
                                chains = chains,
                                WAIC = WAIC,
                                parallel = parallel,
                                cores = cores,
                                seed = seed),
                                bkgd.points = cov.df[samp, c(1:2)])


    if(Distance == "mahalanobis"){
      class(ret.list) <- c("ppmve", Distance)
    }

    if(Distance == "euclidean"){
      class(ret.list) <- c("ppmve", Distance)
    }

    return(ret.list)
  }
}
