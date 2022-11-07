library(mvtnorm)

# Simulation Studies
### sample size
n <- 200
### varying dimensions and different numbers of true regression coefficients.
genr_p_and_active <- function(design){
  # Design 1 (Baseline Design): The number of predictors p=500 and |γ0|=4.
  if (design==1) p <- 500; num_active <- 4
  # Design 2 (Dense model design): The number of predictors p=500 and |γ0|=8.
  if (design==2) p <- 500; num_active <- 8
  # Design 3 (High-dimensional design): The number of predictors p=1000 and |γ0|=4.
  if (design==3) p <- 1000; num_active <- 4
  
  return(c(p, num_active))
}
### Generating X
genr_X <- function(n, design, case){
  p <- genr_p_and_active(design)[1]
  cs_cor <- function(n, rho) matrix(data=rho, nrow=n, ncol=n)+diag(n)*(1-rho)
  ar1_cor <- function(n, rho) rho^(abs(matrix(1:n-1, nrow=n, ncol=n, byrow=TRUE)-(1:n-1)))

  # Case 1: Isotropic design, where Σ = Ip.
  if (case==1) sigm <- diag(p)
  # Case 2: Compound Symmetry Design, where Σij=0.5, if i/=j and Σii=1, for all 1≤i≤j≤p.
  if (case==2) sigm <- cs_cor(p, 0.5)
  # Case 3: Autoregressive Correlated Design; where Σij=0.5|i−j|, for all 1≤i≤j≤p.
  if (case==3) sigm <- ar1_cor(p, 0.5)
  
  set.seed(76)
  X <- rmvnorm(n=n, mean=rep(0, p), sigma=sigm)
  
  return(X)
}
### Generating Beta
genr_Beta <- function(design, setting){
  p <- genr_p_and_active(design)[1]
  num_active <- genr_p_and_active(design)[2]
  
  set.seed(76)
  beta <- rep(0, p)
  # Setting 1: All the entries of β0,γ0 are generated from Unif(0.5, 1.5).
  if (setting==1) beta[1:num_active] <- runif(n=num_active, min=0.5, max=1.5)
  # Setting 2: All the entries of β0, γ0 are set to 1.5.
  if (setting==2) beta[1:num_active] <- 1.5
  # Setting 3: All the entries of β0,γ0 are generated from Unif(1.5, 3).
  if (setting==3) beta[1:num_active] <- runif(n=num_active, min=1.5, max=3)
  # Setting 4: All the entries of β0,γ0 are set to 3.
  if (setting==4) beta[1:num_active] <- 3
  
  return(beta)
}
### Generating Y
# Y <- ifelse(X%*%beta>=0, 1, 0)

################################################################################
### Table 1: Design 1, Case 1 with Setting 1-4
################################################################################
X11 <- genr_X(n, design=1, case=1)

beta11 <- genr_Beta(design=1, setting=1)
Y111 <- ifelse(X11%*%beta11>=0, 1, 0)
df111 <- cbind(X11, Y111)
write.table(df111, file= 'df111.txt', row.names=F, col.names=F)