LocalLocalMahal <- nimble::nimbleCode({
  for(i in 1:n.clim){
    centroid.pres[i] ~ dnorm(cent.mean[i], cent.prec[i])
  }
  
  tau.pres[1:n.clim, 1:n.clim] ~ dwish(R[1:n.clim, 1:n.clim], n.clim + 1)
  
  beta ~ dnorm(beta.mean, beta.prec)
  
  #Centroid likelihood
  for(i in 1:n.data){
    dists.pres[i, 1:n.clim] <- (clim[i, 1:n.clim] - centroid.pres[1:n.clim]) %*% tau.pres[1:n.clim, 1:n.clim] * (clim[i, 1:n.clim] - centroid.pres[1:n.clim]) 
    
    dist.pres[i] <- sum(dists.pres[i, 1:n.clim])
    
    mu[i] <- w[i] * exp((beta - dist.pres[i]/2))
    lambda[i] ~ dpois(mu[i])
  }
})