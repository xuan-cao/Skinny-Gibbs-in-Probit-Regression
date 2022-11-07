rm(list=ls())

source("exact_gibbs_probit.R")
source("skinny_gibbs_probit.R")
source("evaluation.R")


# Import Data
df <- read.table("df111.txt")
X <- scale(data.matrix(df[, -ncol(df)]))
Y <- df[, ncol(df)]
n <- nrow(X)
p <- ncol(X)

# True Gamma 
true_gamma <- c(rep(1, 4), rep(0, 500-4))

# Initializing
Z <- rep(1, n)
tau0_2 <- 1/n
tau1_2 <- p^(2+2*0.0001)/n*10^(-3) #tuning  
q <- p^(-0.1)
gamma <- rep(0, p)
gamma[sample(1:p,20)] = 1

###########################
### Exact Gibbs Sampler ###
###########################
out_exact <- exactGibbs(X, Y, Z, gamma, tau0_2, tau1_2, q)

criteria_exact <- evaluation(actual=true_gamma, predicted=1-(out_exact$outgamma<0.5),
                             X, Y, beta=out_exact$outbeta)
criteria_exact

############################
### Skinny Gibbs Sampler ###
############################
out_skinny <- skinnyGibbs(X, Y, Z, gamma, tau0_2, tau1_2, q)

criteria_skinny <- evaluation(actual=true_gamma, predicted=1-(out_skinny$outgamma<0.5), 
                              X, Y, beta=out_skinny$outbeta)
criteria_skinny