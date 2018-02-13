#################################################################
# Author: M van Smeden                                          #
# Illustrative simulations of small sample bias of log(OR)      #
# Date: Feb 13 2018                                             #
#################################################################

# dependencies #
require(MASS)
require(ggplot2)
require(logistf)
# ----------- #

# input #
P <- 4 # number of covariables (assume the first is exposure, the others are confounders)
N <- 50 # simulated sample size
OR <- rep(2,4) # true odds ratios
corX <- .1 # true (equal) pairwise correlations between covariables
R <- 100 # number of simulation replications, N = 1000 may take a while to run
set.seed(2018)
# ---- #

  sep_warn <- c()
  outputML <- matrix(,ncol=8,nrow=R)
  colnames(outputML) <- c("r","cumN","cumest","iterest","consistency","bias","sqrd_error","MSE")
  outputFR <- outputML
 
  sigma <- matrix(corX,ncol=P,nrow=P) 
    diag(sigma) <- 1
  
  X <- mvrnorm(n=N,mu=rep(0,P),Sigma=sigma)
  Y <- rbinom(N,1,1/(1+exp(-t(log(OR))%*%t(X))))

for(r in 1:R){
  newX <- mvrnorm(n=N,mu=rep(0,P),Sigma=sigma)
  newY <- rbinom(N,1,1/(1+exp(-t(log(OR))%*%t(newX))))
  X <- rbind(X,newX)
  Y <- c(Y,newY)
  
  MLiter <- glm(newY~newX,family="binomial")
  FRiter <- logistf(newY~newX,firth = T) 
   sep_warn[r] <- any(sqrt(diag(vcov(MLiter)))>20)
  MLcum <- glm(Y~X,family="binomial")
  FRcum <- logistf(Y~X,firth = T) 
  
  outputML[r,"r"] <- r ; outputFR[r,"r"] <- r
  outputML[r,"cumN"] <- length(Y);  outputFR[r,"cumN"] <- length(Y)
  outputML[r,"cumest"] <- as.numeric(coef(MLcum)["X1"]); outputFR[r,"cumest"] <- as.numeric(coef(FRcum)["X1"])
  outputML[r,"iterest"] <- as.numeric(coef(MLiter)["newX1"]); outputFR[r,"iterest"] <- as.numeric(coef(FRiter)["newX1"])
  outputML[r,"consistency"] <- outputML[r,"cumest"]-log(OR[1]); outputFR[r,"consistency"] <- outputFR[r,"cumest"]-log(OR[1])
  outputML[r,"bias"] <- mean(outputML[,"iterest"],na.rm=T)-log(OR[1]); outputFR[r,"bias"] <-  mean(outputFR[,"iterest"],na.rm=T)-log(OR[1])
  outputML[r,"sqrd_error"] <- (outputML[r,"iterest"]-log(OR[1]))^2; outputFR[r,"sqrd_error"] <- (outputFR[r,"iterest"]-log(OR[1]))^2
  outputML[r,"MSE"] <- mean(outputML[,"sqrd_error"],na.rm=T); outputFR[r,"MSE"] <- mean(outputFR[,"sqrd_error"],na.rm=T)
  
  print(paste(r,"of",R))
}

table(sep_warn)  # note: if this table indicates all "FALSE", separation might have occured. Simulation results might then not be trustworthy.

  outputFR[R,"bias"]/log(OR[1])
  outputML[R,"bias"]/log(OR[1])
  