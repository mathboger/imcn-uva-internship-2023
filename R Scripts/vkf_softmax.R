# Script to test the VKF + softmax model (referred as stepping stone model in the report)
# Mathematical definitions and parameter recovery (data & plots)

rm(list=ls())
library(rtdists)
library(DEoptim)
library(ggplot2)
library(gridExtra)

### Generate data following the predictive model of the binary VKF
### Connect to Softmax with M (M*beta) (i.e., M serves as a Q-value)

# VKF parameters
m0 <- 0.5 # Initial state
v0 <- 1 # Initial volatility
lambda <- 0.5 # Volatility update rate
omega <- 1 # Observation noise

# Softmax parameter
beta <- 3 # Inverse temperature parameter

# Experiment parameters
nTrials <- 100
nChoices <- 2
reversalPoint <- 70
pReward <- c(.7, .3)
nSubjects <- 1 # Parameter Recovery
nIterations <- 50 # Parameter Recovery
# nSubjects is nested inside nIterations, with the parameters for parameter recovery
# being randomly sampled and renewed every iteration and kept constant for the number
# of participants inside one iteration.

# Function that takes in VKF, softmax, and experiment parameters and outputs simulated data along the model
generateData <- function(nTrials, nChoices, reversalPoint, pReward, m0, v0, lambda, omega, beta) {
  
  # Initializing variables
  # VKF
  K <- matrix(NA, nrow=nTrials, ncol=nChoices)
  Alpha <- matrix(NA, nrow=nTrials, ncol=nChoices)
  M <- matrix(NA, nrow=nTrials, ncol=nChoices)
  W <- matrix(NA, nrow=nTrials, ncol=nChoices)
  CoW <- matrix(NA, nrow=nTrials, ncol=nChoices)
  V <- matrix(NA, nrow=nTrials, ncol=nChoices)
  # Output
  choice <- numeric(nTrials)
  reward <- logical(nTrials)
  
  K[1,] <- rep(v0/(v0 + omega), nChoices)
  Alpha[1,] <- rep(sqrt(v0), nChoices)
  M[1,] <- rep(m0, nChoices)
  W[1,] <- (1 - K[1,])*v0
  CoW[1,] <- rep(0, nChoices)
  V[1,] <- v0 + lambda*(W[1,] - v0)
  
  # Generate data
  for(i in 2:nTrials) {
    
    # Check for reversal
    if(i == reversalPoint) {
      pReward <- c(pReward[2], pReward[1])
    }
    
    # Make observation
    pChoice <- exp(beta*M[i-1,])/sum(exp(beta*M[i-1,])) # Softmax function
    choice[i] <- rbinom(1, 1, pChoice[1]) # Draw a choice
    choice[i] <- ifelse(choice[i]==1, 1, 2) # Recode from 0/1 to 1/2
    reward[i] <- runif(1) < pReward[choice[i]] # Draw a reward
    
    # Update variables (only change the values of the VKF associated to the choice made)
    K[i,] <- K[i-1,]
    K[i, choice[i]] <- (W[i-1, choice[i]] + V[i-1, choice[i]])/(W[i-1, choice[i]] + V[i-1, choice[i]] + omega)
    Alpha[i,] <- Alpha[i-1,]
    Alpha[i, choice[i]] <- sqrt(W[i-1, choice[i]] + V[i-1, choice[i]])
    M[i,] <- M[i-1,]
    M[i, choice[i]] <- M[i-1, choice[i]] + Alpha[i, choice[i]]*(reward[i] - 1/(1 + exp(-M[i-1, choice[i]])))
    W[i,] <- W[i-1,]
    W[i, choice[i]] <- (1 - K[i, choice[i]])*(W[i-1, choice[i]] + V[i-1, choice[i]])
    CoW[i,] <- CoW[i-1,]
    CoW[i, choice[i]] <- (1 - K[i, choice[i]])*W[i-1, choice[i]]
    V[i,] <- V[i-1,]
    V[i, choice[i]] <- V[i-1, choice[i]] + lambda*((M[i, choice[i]] - M[i-1, choice[i]])**2 + W[i-1, choice[i]] + W[i, choice[i]] - 2*CoW[i, choice[i]] - V[i-1, choice[i]])
  }
  
  model <- data.frame(K=K, Alpha=Alpha, M=M, W=W, CoW=CoW, V=V, choice=choice, reward=reward)
  return(model)
}

# Usage example
model <- generateData(nTrials, nChoices, reversalPoint, pReward, m0, v0, lambda, omega, beta)

### Likelihood

# Function that takes in the observed outcome and the proposed model parameters and outputs
# the likelihood of that model having generated the data or the model itself (intermediate variables)
ll_func <- function(choice, reward, nTrials, nChoices, m0, v0, lambda, omega, beta, returnModel=FALSE) {
  
  # Initialize variables
  K <- matrix(NA, nrow=nTrials, ncol=nChoices)
  Alpha <- matrix(NA, nrow=nTrials, ncol=nChoices)
  M <- matrix(NA, nrow=nTrials, ncol=nChoices)
  W <- matrix(NA, nrow=nTrials, ncol=nChoices)
  CoW <- matrix(NA, nrow=nTrials, ncol=nChoices)
  V <- matrix(NA, nrow=nTrials, ncol=nChoices)
  
  K[1,] <- rep(v0/(v0 + omega), nChoices)
  Alpha[1,] <- rep(sqrt(v0), nChoices)
  M[1,] <- rep(m0, nChoices)
  W[1,] <- (1 - K[1,])*v0
  CoW[1,] <- rep(0, nChoices)
  V[1,] <- v0 + lambda*(W[1,] - v0)
  
  for(i in 2:nTrials) {
    
    # Update variables (only change the values of the VKF associated to the choice made)
    K[i,] <- K[i-1,]
    K[i, choice[i]] <- (W[i-1, choice[i]] + V[i-1, choice[i]])/(W[i-1, choice[i]] + V[i-1, choice[i]] + omega)
    Alpha[i,] <- Alpha[i-1,]
    Alpha[i, choice[i]] <- sqrt(W[i-1, choice[i]] + V[i-1, choice[i]])
    M[i,] <- M[i-1,]
    M[i, choice[i]] <- M[i-1, choice[i]] + Alpha[i, choice[i]]*(reward[i] - 1/(1 + exp(-M[i-1, choice[i]])))
    W[i,] <- W[i-1,]
    W[i, choice[i]] <- (1 - K[i, choice[i]])*(W[i-1, choice[i]] + V[i-1, choice[i]])
    CoW[i,] <- CoW[i-1,]
    CoW[i, choice[i]] <- (1 - K[i, choice[i]])*W[i-1, choice[i]]
    V[i,] <- V[i-1,]
    V[i, choice[i]] <- V[i-1, choice[i]] + lambda*((M[i, choice[i]] - M[i-1, choice[i]])**2 + W[i-1, choice[i]] + W[i, choice[i]] - 2*CoW[i, choice[i]] - V[i-1, choice[i]])
  }
  
  # Calculate likelihood
  LL <- numeric(nTrials)
  for (i in 2:nTrials) {
    pChoice <- exp(beta*M[i-1,])/sum(exp(beta*M[i-1,]))
    LL[i] <- pChoice[choice[i]]
  }
  LL[1] <- 1 # Placeholder that won't change the value of sum(log()) but will make LL the right size to couple it to a model
  LL[LL == 0] <- 1e-20 # Adjustment for log in case we get a misfit
  
  if(returnModel) {
    model <- data.frame(K=K, Alpha=Alpha, M=M, W=W, CoW=CoW, V=V, choice=choice, reward=reward)
    return(model)
  }
  return(LL)
}

### Profile Data Generation

rangeM0 <- seq(0.01, 5, 0.01)
rangeV0 <- seq(0.01, 5, 0.01)
rangeLambda <- seq(0.01, 0.99, 0.01)
rangeOmega <- seq(0.01, 5, 0.01)
rangeBeta <- seq(0.01, 5, 0.01)

LLM0 <- numeric(length(rangeM0))
for(i in 1:length(rangeM0)) {
  LLM0[i] <- sum(log(ll_func(model$choice, model$reward, nTrials, nChoices, rangeM0[i], v0, lambda, omega, beta)))
}
LLV0 <- numeric(length(rangeV0))
for(i in 1:length(rangeV0)) {
  LLV0[i] <- sum(log(ll_func(model$choice, model$reward, nTrials, nChoices, m0, rangeV0[i], lambda, omega, beta)))
}
LLLambda <- numeric(length(rangeLambda))
for(i in 1:length(rangeLambda)) {
  LLLambda[i] <- sum(log(ll_func(model$choice, model$reward, nTrials, nChoices, m0, v0, rangeLambda[i], omega, beta)))
}
LLOmega <- numeric(length(rangeOmega))
for(i in 1:length(rangeOmega)) {
  LLOmega[i] <- sum(log(ll_func(model$choice, model$reward, nTrials, nChoices, m0, v0, lambda, rangeOmega[i], beta)))
}
LLBeta <- numeric(length(rangeBeta))
for(i in 1:length(rangeBeta)) {
  LLBeta[i] <- sum(log(ll_func(model$choice, model$reward, nTrials, nChoices, m0, v0, lambda, omega, rangeBeta[i])))
}
save(rangeM0, rangeV0, rangeLambda, rangeOmega, rangeBeta, LLM0, LLV0, LLLambda, LLOmega, LLBeta, m0, v0, lambda, omega, beta, file="../Data/Simulated Data/Stepping Stone Model Profile/profile.RData")

### Parameter Recovery Data Generation

# Wrapper function for DEoptim
wrapper <- function(pars, nTrials, nChoices, trueModel) {
  LL <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, pars[1], pars[2], pars[3], pars[4], pars[5])
  return(-sum(log(LL)))
}

startTime <- Sys.time() # Keep track of how long it takes
for(i in 1:nIterations) {
  # Get save file / directory
  save_dir = paste0("../Data/Simulated Data/Stepping Stone Model Parameter Recovery/iteration-", i)
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive=TRUE)
  
  # Sample new ground truth
  trueM0 <- runif(1, 0.01, 5)
  trueV0 <- runif(1, 0.01, 5)
  trueLambda <- runif(1, 0.01, 0.99)
  trueOmega <- runif(1, 0.01, 5)
  trueBeta <- runif(1, 0.01, 5)
  
  for(j in 1:nSubjects) {
    
    # Get true model and likelihood
    trueModel <- generateData(nTrials, nChoices, reversalPoint, pReward, trueM0, trueV0, trueLambda, trueOmega, trueBeta)
    trueLL <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, trueM0, trueV0, trueLambda, trueOmega, trueBeta)
    
    # Get fit model and likelihood
    out <- DEoptim(wrapper, lower=c(0.01, 0.01, 0.01, 0.01, 0.01), upper=c(5, 5, 0.99, 5, 5), nTrials=nTrials, nChoices=nChoices, trueModel=trueModel)
    fitM0 <- out$optim$bestmem[1]
    fitV0 <- out$optim$bestmem[2]
    fitLambda <- out$optim$bestmem[3]
    fitOmega <- out$optim$bestmem[4]
    fitBeta <- out$optim$bestmem[5]
    fitModel <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, fitM0, fitV0, fitLambda, fitOmega, fitBeta, returnModel=TRUE)
    fitLL <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, fitM0, fitV0, fitLambda, fitOmega, fitBeta)
    
    # Save ground truth and fit to save file
    save(trueM0, trueV0, trueLambda, trueOmega, trueBeta, 
         fitM0, fitV0, fitLambda, fitOmega, fitBeta,
         trueModel, fitModel, trueLL, fitLL, file=paste0(save_dir, "/subject-", j, ".RData"))
  }
}
endTime <- Sys.time()
print(endTime-startTime) # Check how long it took

### Plotting

# Figure 3 of the report
# Fetch data
load("../Data/Simulated Data/Stepping Stone Model Profile/profile.RData")

# Plot
m0PP <- ggplot(data.frame(rangeM0=rangeM0, LLM0=LLM0), aes(x=rangeM0, y=LLM0)) +
  theme_classic() + ggtitle(expression((A)~m[0])) +
  geom_line() + xlab(expression(m[0])) + ylab("Likelihood") +
  geom_vline(xintercept=m0, col="blue")

v0PP <- ggplot(data.frame(rangeV0=rangeV0, LLV0=LLV0), aes(x=rangeV0, y=LLV0)) +
  theme_classic() + ggtitle(expression((B)~v[0])) +
  geom_line() + xlab(expression(v[0])) + ylab("Likelihood") +
  geom_vline(xintercept=v0, col="blue")

lambdaPP <- ggplot(data.frame(rangeLambda=rangeLambda, LLLambda=LLLambda), aes(x=rangeLambda, y=LLLambda)) +
  theme_classic() + ggtitle("(C) λ") +
  geom_line() + xlab("λ") + ylab("Likelihood") +
  geom_vline(xintercept=lambda, col="blue")

omegaPP <- ggplot(data.frame(rangeOmega=rangeOmega, LLOmega=LLOmega), aes(x=rangeOmega, y=LLOmega)) +
  theme_classic() + ggtitle("(D) ω") +
  geom_line() + xlab("ω") + ylab("Likelihood") +
  geom_vline(xintercept=omega, col="blue")

betaPP <- ggplot(data.frame(rangeBeta=rangeBeta, LLBeta=LLBeta), aes(x=rangeBeta, y=LLBeta)) +
  theme_classic() + ggtitle("(E) β") +
  geom_line() + xlab("β") + ylab("Likelihood") +
  geom_vline(xintercept=beta, col="blue")

fig3 <- arrangeGrob(m0PP, v0PP, lambdaPP, omegaPP, betaPP, ncol=2, nrow=3, top="Stepping Stone Model Profile Plot")
ggsave("../Figures/fig3.png", plot=fig3, width = 2 * 400, height = 3 * 400 + 20, units = "px", dpi=72)

# Figure 6 of the report
# Fetch data
trueM0s <- numeric(nIterations*nSubjects)
trueV0s <- numeric(nIterations*nSubjects)
trueLambdas <- numeric(nIterations*nSubjects)
trueOmegas <- numeric(nIterations*nSubjects)
trueBetas <- numeric(nIterations*nSubjects)
trueLLs <- numeric(nIterations*nSubjects)
fitM0s <- numeric(nIterations*nSubjects)
fitV0s <- numeric(nIterations*nSubjects)
fitLambdas <- numeric(nIterations*nSubjects)
fitOmegas <- numeric(nIterations*nSubjects)
fitBetas <- numeric(nIterations*nSubjects)
fitLLs <- numeric(nIterations*nSubjects)
k <- 1
for(i in 1:nIterations) {
  save_dir <- paste0("../Data/Simulated Data/Stepping Stone Model Parameter Recovery/iteration-", i)
  for(j in 1:nSubjects) {
    # Load data
    load(paste0(save_dir, "/subject-", j, ".RData"))
    trueM0s[k] <- trueM0
    trueV0s[k] <- trueV0
    trueLambdas[k] <- trueLambda
    trueOmegas[k] <- trueOmega
    trueBetas[k] <- trueBeta
    trueLLs[k] <- -sum(log(trueLL))
    fitM0s[k] <- fitM0
    fitV0s[k] <- fitV0
    fitLambdas[k] <- fitLambda
    fitOmegas[k] <- fitOmega
    fitBetas[k] <- fitBeta
    fitLLs[k] <- -sum(log(fitLL))
    k <- k + 1
  }
}

# Plot
lmM0 <- lm(fitM0s ~ trueM0s)
corM0 <- cor(trueM0s, fitM0s, method="spearman")
m0Plot <- ggplot(data.frame(trueM0s=trueM0s, fitM0s=fitM0s), aes(x=trueM0s, y=fitM0s)) +
  theme_classic() + ggtitle(expression((A)~m[0]~recovery)) +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab(expression(True~m[0])) + ylab(expression(Estimated~m[0])) +
  geom_abline(aes(intercept=lmM0$coefficients[1], slope=lmM0$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") + 
  geom_text(x=5, y=0, label=paste0("r = ", round(corM0, 2)))

lmV0 <- lm(fitV0s ~ trueV0s)
corV0 <- cor(trueV0s, fitV0s, method="spearman")
v0Plot <- ggplot(data.frame(trueV0s=trueV0s, fitV0s=fitV0s), aes(x=trueV0s, y=fitV0s)) +
  theme_classic() + ggtitle(expression((B)~v[0]~recovery)) +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab(expression(True~v[0])) + ylab(expression(Estimated~v[0])) +
  geom_abline(aes(intercept=lmV0$coefficients[1], slope=lmV0$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=5, y=0, label=paste0("r = ", round(corV0, 2)))

lmLambda <- lm(fitLambdas ~ trueLambdas)
corLambda <- cor(trueLambdas, fitLambdas, method="spearman")
lambdaPlot <- ggplot(data.frame(trueLambdas=trueLambdas, fitLambdas=fitLambdas), aes(x=trueLambdas, y=fitLambdas)) +
  theme_classic() + ggtitle("(C) λ recovery") +
  geom_point() + xlim(0, 1) + ylim(0, 1) + xlab("True λ") + ylab("Estimated λ") +
  geom_abline(aes(intercept=lmLambda$coefficients[1], slope=lmLambda$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=1, y=0, label=paste0("r = ", round(corLambda, 2)))

lmOmega <- lm(fitOmegas ~ trueOmegas)
corOmega <- cor(trueOmegas, fitOmegas, method="spearman")
omegaPlot <- ggplot(data.frame(trueOmegas=trueOmegas, fitOmegas=fitOmegas), aes(x=trueOmegas, y=fitOmegas)) +
  theme_classic() + ggtitle("(D) ω recovery") +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab("True ω") + ylab("Estimated ω") +
  geom_abline(aes(intercept=lmOmega$coefficients[1], slope=lmOmega$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=5, y=0, label=paste0("r = ", round(corOmega, 2)))

lmBeta <- lm(fitBetas ~ trueBetas)
corBeta <- cor(trueBetas, fitBetas, method="spearman")
betaPlot <- ggplot(data.frame(trueBetas=trueBetas, fitBetas=fitBetas), aes(x=trueBetas, y=fitBetas)) +
  theme_classic() + ggtitle("(E) β recovery") +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab("True β") + ylab("Estimated β") +
  geom_abline(aes(intercept=lmBeta$coefficients[1], slope=lmBeta$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=5, y=0, label=paste0("r = ", round(corBeta, 2)))

lmLL <- lm(fitLLs ~ trueLLs)
corLL <- cor(trueLLs, fitLLs, method="spearman")
LLPlot <- ggplot(data.frame(trueLLs=trueLLs, fitLLs=fitLLs), aes(x=trueLLs, y=fitLLs)) +
  theme_classic() + ggtitle("(F) Likelihood recovery") +
  geom_point() + xlim(0, 70) + ylim(0, 70) + xlab("True likelihood") + ylab("Estimated likelihood") +
  geom_abline(aes(intercept=lmLL$coefficients[1], slope=lmLL$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=70, y=0, label=paste0("r = ", round(corLL, 2)))

fig6 <- arrangeGrob(m0Plot, v0Plot, lambdaPlot, omegaPlot, betaPlot, LLPlot, ncol=2, nrow=3, top="Stepping Stone Model Parameter Recovery")
ggsave("../Figures/fig6.png", plot=fig6, width = 2 * 400, height = 3 * 400 + 20, units = "px", dpi=72)