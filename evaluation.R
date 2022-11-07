#rm(list=ls())
#setwd("/Users/jrouyang/Desktop/skinny_gibbs_probit/oyang/10_03M")


################################################################################
### Define evaluation Function
################################################################################

evaluation <- function(actual, predicted, beta, X, Y){
  
  true.idx <- which(actual==1)
  false.idx <- which(actual==0)
  positive.idx <- which(predicted==1)
  negative.idx <- which(predicted==0)
  
  TP <- length(intersect(true.idx, positive.idx))
  FP <- length(intersect(false.idx, positive.idx))
  FN <- length(intersect(true.idx, negative.idx))
  TN <- length(intersect(false.idx, negative.idx))
  
  
  Sensitivity <- TP/(TP+FN)
  if ((TP+FN)==0) Sensitivity <- 1

  Specific <- TN/(TN+FP)
  if ((TN+FP)==0) Specific <- 1
  
  MCC.denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if ((TN+FP)==0) MCC <- 1
  
  MSPE <- mean((pnorm(X%*%beta)-Y)^2)
  
  
  return(list(Sensitivity=Sensitivity, Specific=Specific, MCC=MCC, MSPE=MSPE, 
              TP=TP, FP=FP, TN=TN, FN=FN))
}


