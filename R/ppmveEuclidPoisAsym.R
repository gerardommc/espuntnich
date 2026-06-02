EuclidAsym <- nimble::nimbleCode({
  
  # Priors for Niche Centroids, Precisions, and Skewness
  for(i in 1:n.clim){
    centroid.pres[i] ~ dnorm(cent.mean[i], cent.prec[i])
    tau.pres[i]      ~ dunif(tau.min[i], tau.max[i])
    alpha[i]         ~ dnorm(alpha.mean[i], alpha.prec[i])
  }
  
  # Prior for Global Intercept
  beta ~ dnorm(beta.mean, beta.tau)
  
  # Asymmetric Centroid Likelihood
  for(i in 1:n.data){
    
    dists.pres[i, 1:n.clim] <- (clim[i, 1:n.clim] - centroid.pres[1:n.clim])^2 * tau.pres[1:n.clim] 
    dist.pres[i]            <- sum(dists.pres[i, 1:n.clim])
    
    skew.dims[i, 1:n.clim]  <- (clim[i, 1:n.clim] - centroid.pres[1:n.clim]) * alpha[1:n.clim]
    skew.val[i]             <- sum(skew.dims[i, 1:n.clim])
    
    mu[i]      <- w[i] * exp(beta - dist.pres[i]/2) * 2 * phi(skew.val[i])
    lambda[i]  ~ dpois(mu[i])
  }
})