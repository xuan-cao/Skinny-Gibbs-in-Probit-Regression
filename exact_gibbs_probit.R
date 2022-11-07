library(mvtnorm)
library(truncnorm)

################################################################################
### Define Exact Gibbs Function
################################################################################
exactGibbs <- function(X, Y, Z, gamma, tau0_2, tau1_2, q, nburn=2000, niter=2000){
  #####
  # X is a nxp matrix
  # Y is a 1xn vector with 0's or 1's
  # Z is a 1xn vector
  # gamma is a 1xp vector with 0's or 1's
  # tau0 is a real number 
  # tau1 is a real number
  # q is a real number
  #####
  n <- nrow(X)         
  p <- ncol(X)
  
  outbeta <- rep(0, p)
  outgamma <- rep(0, p)
  outZ <- rep(0, n)
  
  inv_tau0_2 <- 1/tau0_2
  inv_tau1_2 <- 1/tau1_2
  const <- (q*sqrt(tau0_2))/((1-q)*sqrt(tau1_2))
  
  # start updating
  for (itr in 1:(nburn+niter)){
    #print(paste("iteration", itr))
    # Check Processing
    if (itr %% 100 == 0) print(paste("finish iteration", itr))
    
    ##### update beta
    D_gamma <- diag(gamma*inv_tau1_2+(1-gamma)*inv_tau0_2)
    sigma1 <- solve(t(X)%*%X+D_gamma)
    mean1 <- sigma1%*%t(X)%*%Y
    R <- chol(sigma1)
    temp <- t(R)%*%matrix(rnorm(1000*p), p)
    beta <- mean1+diag(temp%*%t(temp)/1000)
    
    ##### update gamma
    for (j in 1:p){
      temp1 <- beta[j]^2*0.5*(inv_tau0_2-inv_tau1_2)
      log_d_j <- log(const)+temp1
      gamma[j] <- rbinom(1, 1, exp(log_d_j)/(1+exp(log_d_j)))
    }

    #### update Z
    for (i in 1:n){
      if (Y[i]==1) Z[i] <- rtruncnorm(1, a=0, mean=X[i,]%*%beta, sd=1)
      if (Y[i]==0) Z[i] <- rtruncnorm(1, b=0, mean=X[i,]%*%beta, sd=1)
    }
    
    if (itr > nburn) {
      outbeta <- outbeta+beta
      outgamma <- outgamma+gamma
      outZ <- outZ+Z
    }
  }
  
  return(list(outbeta=outbeta/niter, outgamma=outgamma/niter, outZ=outZ/niter))
}


