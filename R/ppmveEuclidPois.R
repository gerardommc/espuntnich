modelCode <- nimbleCode({
  for(i in 1:n.clim){
    centroid.pres[i] ~ dnorm(0, 0.1)
    tau.pres[i]  ~ dunif(0, 100)
  }
  
  beta ~ dnorm(0, 1.0E-3)
  
  #Centroid likelihood
  for(i in 1:n.data){
    dists.pres[i, 1:n.clim] <- (clim[i, 1:n.clim] - centroid.pres[1:n.clim])^2 * tau.pres[1:n.clim] 
    
    dist.pres[i] <- sum(dists.pres[i, 1:n.clim])
    
    mu[i] <-  w[i] * exp(beta - dist.pres[i]/2)
    lambda[i] ~ dpois(mu[i])
  }
})