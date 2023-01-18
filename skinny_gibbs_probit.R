library(mvtnorm)
library(truncnorm)

################################################################################
### Define skinny Gibbs Function
################################################################################
skinnyGibbs <- function(X, Y, Z, gamma, tau0_2, tau1_2, q, nburn=2000, niter=2000){
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
    idx_active <- which(gamma==1)
    idx_inactive <- which(gamma==0)
    n_active <- length(idx_active)
    
    beta <- rep(0, p)
    if (n_active==0) beta <- rnorm(n=p, mean=0, sd=1/sqrt(n-1+inv_tau0_2)) 
    else {
      # Generate active beta
      X_active <- matrix(X[, idx_active], nrow=n)
      sigma1 <- solve(t(X_active)%*%X_active+inv_tau1_2*diag(n_active))
      mean1 <- sigma1%*%t(X_active)%*%Z
      beta[idx_active] <- rmvnorm(n=1, mean=mean1, sigma=sigma1)
      # Generate inactive beta 
      beta[idx_inactive] <- rnorm(n=p-n_active, mean=0, sd=1/sqrt(n-1+inv_tau0_2))
    }      
    
    ##### update gamma
    if (n_active==0) temp2_2_2 <- Z
    else{
      X_active <- matrix(X[, idx_active], nrow=n) # guarantee it is matrix form. 
      beta_active <- beta[idx_active]
      temp2_2_2 <- Z-X_active%*%beta_active # dim(temp2_2_2) is nx1
    }
    for (j in 1:p){
      temp1 <- beta[j]^2*0.5*(inv_tau0_2-inv_tau1_2)
      temp2_1 <- beta[j]*X[,j]
      if (j %in% idx_active){
        temp2_2_1 <- Z-matrix(X_active[, -which(idx_active==j)], 
                              nrow=n)%*%beta_active[-which(idx_active==j)]
        temp2 <- temp2_1%*%temp2_2_1
      }
      else temp2 <- temp2_1%*%temp2_2_2
      
      log_d_j_star <- log(const)+temp1+temp2
      if (exp(log_d_j_star)==Inf) gamma[j] <- 1
      else gamma[j] <- rbinom(1, 1, exp(log_d_j_star)/(1+exp(log_d_j_star)))
    }
    
    #### update Z
    idx_active <- which(gamma==1)
    X_active <- matrix(X[, idx_active], nrow=n) # guarantee it is matrix form. 
    beta_active <- beta[idx_active]
    for (i in 1:n){
      if (Y[i]==1) Z[i] <- rtruncnorm(1, a=0, mean=X_active[i,]%*%beta_active, sd=1)
      if (Y[i]==0) Z[i] <- rtruncnorm(1, b=0, mean=X_active[i,]%*%beta_active, sd=1)
    }
    
    if (itr > nburn) {
      outbeta <- outbeta+beta
      outgamma <- outgamma+gamma
      outZ <- outZ+Z
    }
  }
  
  return(list(outbeta=outbeta/niter, outgamma=outgamma/niter, outZ=outZ/niter))
}
