GlobalMahal <- nimble::nimbleCode({
  for(i in 1:n.clim){
    centroid.pres[i] ~ dnorm(0, 0.1)
    mu.back[i]  ~ dnorm(0, 1.0E-4)
  }
  
  tau.pres[1:n.clim, 1:n.clim] ~ dwish(R[1:n.clim, 1:n.clim], n.clim + 1)
  
  for(j in 1:n.back){
    clim[j, 1:n.clim] ~ dmnorm(mu.back[1:n.clim], tau.pres[1:n.clim, 1:n.clim])
  }
  
  beta ~ dnorm(0, 1.0E-3)
  
  #Centroid likelihood
  for(i in 1:n.data){
    dists.pres[i, 1:n.clim] <- (clim[i, 1:n.clim] - centroid.pres[1:n.clim]) %*% tau.pres[1:n.clim, 1:n.clim] * (clim[i, 1:n.clim] - centroid.pres[1:n.clim]) 
    
    dist.pres[i] <- sum(dists.pres[i, 1:n.clim])
    
    mu[i] <- w[i] * exp((beta - dist.pres[i]/2))
    lambda[i] ~ dpois(mu[i])
  }
})