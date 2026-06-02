LocalLocalMahalAsym <- nimble::nimbleCode({
  for(i in 1:n.clim){
    centroid.pres[i] ~ dnorm(cent.mean[i], cent.prec[i])
    alpha[i] ~ dnorm(alpha.mean[i], alpha.prec[i]) # Asymmetric shape parameter prior
  }
  
  tau.pres[1:n.clim, 1:n.clim] ~ dwish(R[1:n.clim, 1:n.clim], n.clim + 1)
  
  beta ~ dnorm(beta.mean, beta.prec)
  
  # Centroid likelihood
  for(i in 1:n.data){
    # Symmetric Mahalanobis distance component
    dists.pres[i, 1:n.clim] <- (clim[i, 1:n.clim] - centroid.pres[1:n.clim]) %*% tau.pres[1:n.clim, 1:n.clim] * (clim[i, 1:n.clim] - centroid.pres[1:n.clim]) 
    dist.pres[i] <- sum(dists.pres[i, 1:n.clim])
    
    # Skewness modulation component
    skew.val[i] <- sum((clim[i, 1:n.clim] - centroid.pres[1:n.clim]) * alpha[1:n.clim])
    
    mu[i] <- w[i] * exp((beta - dist.pres[i]/2)) * 2 * phi(skew.val[i])
    lambda[i] ~ dpois(mu[i])
  }
})