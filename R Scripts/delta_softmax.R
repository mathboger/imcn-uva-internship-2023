# Script to test the simple delta rule + softmax model (referred as baseline model in the report)
# Mathematical definitions and parameter recovery (data & plots)

rm(list=ls())
library(rtdists)
library(DEoptim)
library(ggplot2)
library(gridExtra)

### Generate data

# Simple delta rule parameter
alpha <- 0.3 # Learning rate

# Softmax parameter
beta <- 3 # Inverse temperature parameter

# Experiment parameters (for data simulation)
nTrials <- 100
nChoices <- 2
reversalPoint <- 70
pReward <- c(.7, .3)
nSubjects <- 1 # Parameter Recovery
nIterations <- 50 # Parameter Recovery
# nSubjects is nested inside nIterations, with the parameters for parameter recovery
# being randomly sampled and renewed every iteration and kept constant for the number
# of participants inside one iteration.

# Function that takes in simple delta rule, softmax, and experiment parameters and outputs simulated data along the model
generateData <- function(nTrials, nChoices, reversalPoint, pReward, alpha, beta) {
  
  # Initializing variables
  Q <- matrix(NA, nrow=nTrials, ncol=nChoices) # Q-value for RL
  # Output
  choice <- numeric(nTrials)
  reward <- logical(nTrials)
  
  Q[1,] <- rep(0, nChoices)
  
  # Generate data
  for(i in 2:nTrials) {
    
    # Check for reversal
    if(i == reversalPoint) {
      pReward <- c(pReward[2], pReward[1])
    }
    
    # Make observation
    pChoice <- exp(beta*Q[i-1,])/sum(exp(beta*Q[i-1,])) # Softmax function
    choice[i] <- rbinom(1, 1, pChoice[1]) # Draw a choice
    choice[i] <- ifelse(choice[i]==1, 1, 2) # Recode from 0/1 to 1/2
    reward[i] <- runif(1) < pReward[choice[i]] # Draw a reward
    
    # Update variables (only change the Q-value associated to the choice made)
    Q[i,] <- Q[i-1,]
    Q[i, choice[i]] <- Q[i-1, choice[i]] + alpha*(reward[i] - Q[i-1, choice[i]])
  }
  
  model <- data.frame(Q=Q, choice=choice, reward=reward)
  return(model)
}

# Usage example
model <- generateData(nTrials, nChoices, reversalPoint, pReward, alpha, beta)

### Likelihood

# Function that takes in the observed outcome and the proposed model parameters and outputs
# the likelihood of that model having generated the data or the model itself (intermediate variables)
ll_func <- function(choice, reward, nTrials, nChoices, alpha, beta, returnModel=FALSE) {
  
  # Initialize variables
  Q <- matrix(NA, nrow=nTrials, ncol=nChoices)
  
  Q[1,] <- rep(0, nChoices)
  
  for(i in 2:nTrials) {
    
    # Update variables (only change the Q-value associated to the choice made)
    Q[i,] <- Q[i-1,]
    Q[i, choice[i]] <- Q[i-1, choice[i]] + alpha*(reward[i] - Q[i-1, choice[i]])
  }
  
  ## Calculate likelihood
  LL <- numeric(nTrials)
  for(i in 2:nTrials) {
    pChoice <- exp(beta*Q[i-1,])/sum(exp(beta*Q[i-1,]))
    LL[i] <- pChoice[choice[i]]
  }
  LL[1] <- 1 # Placeholder that won't change the value of sum(log()) but will make LL the right size to couple it to a model
  LL[LL == 0] <- 1e-20 # Adjustment for log in case we get a misfit
  
  if(returnModel) {
    model <- data.frame(Q=Q, choice=choice, reward=reward)
    return(model)
  }
  return(LL)
}

### Profile Data Generation

rangeAlpha <- seq(0.01, 1, 0.01)
rangeBeta <- seq(0.01, 5, 0.01)

LLAlpha <- numeric(length(rangeAlpha))
for(i in 1:length(rangeAlpha)) {
  LLAlpha[i] <- sum(log(ll_func(model$choice, model$reward, nTrials, nChoices, rangeAlpha[i], beta)))
}
LLBeta <- numeric(length(rangeBeta))
for(i in 1:length(rangeBeta)) {
  LLBeta[i] <- sum(log(ll_func(model$choice, model$reward, nTrials, nChoices, alpha, rangeBeta[i])))
}
save(rangeAlpha, rangeBeta, LLAlpha, LLBeta, alpha, beta, file="../Data/Simulated Data/Baseline Model Profile/profile.RData")

### Parameter Recovery Data Generation

# Wrapper function for DEoptim
wrapper <- function(pars, nTrials, nChoices, trueModel) {
  LL <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, pars[1], pars[2])
  return(-sum(log(LL)))
}

startTime <- Sys.time() # Keep track of how long it takes
for(i in 1:nIterations) {
  # Get save file / directory
  save_dir = paste0("../Data/Simulated Data/Baseline Model Parameter Recovery/iteration-", i)
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive=TRUE)
  
  # Sample new ground truth
  trueAlpha <- runif(1, 0.01, 1)
  trueBeta <- runif(1, 0.01, 5)
  
  for(j in 1:nSubjects) {
    
    # Get true model and likelihood
    trueModel <- generateData(nTrials, nChoices, reversalPoint, pReward, trueAlpha, trueBeta)
    trueLL <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, trueAlpha, trueBeta)
    
    # Get fit model and likelihood
    out <- DEoptim(wrapper, lower=c(0.01, 0.01), upper=c(1, 5), nTrials=nTrials, nChoices=nChoices, trueModel=trueModel)
    fitAlpha <- out$optim$bestmem[1]
    fitBeta <- out$optim$bestmem[2]
    fitModel <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, fitAlpha, fitBeta, returnModel=TRUE)
    fitLL <- ll_func(trueModel$choice, trueModel$reward, nTrials, nChoices, fitAlpha, fitBeta)
    
    # Save ground truth and fit to save file
    save(trueAlpha, trueBeta,
         fitAlpha, fitBeta,
         trueModel, fitModel, trueLL, fitLL, file=paste0(save_dir, "/subject-", j, ".RData"))
  }
}
endTime <- Sys.time()
print(endTime-startTime) # Check how long it took

### Plotting

# Figure 2 of the report
# Fetch data
load("../Data/Simulated Data/Baseline Model Profile/profile.RData")

# Plot
alphaPP <- ggplot(data.frame(rangeAlpha=rangeAlpha, LLAlpha=LLAlpha), aes(x=rangeAlpha, y=LLAlpha)) +
  theme_classic() + ggtitle("(A) α") +
  geom_line() + xlab("α") + ylab("Likelihood") +
  geom_vline(xintercept=alpha, col="blue")

betaPP <- ggplot(data.frame(rangeBeta=rangeBeta, LLBeta=LLBeta), aes(x=rangeBeta, y=LLBeta)) +
  theme_classic() + ggtitle("(B) β") +
  geom_line() + xlab("β") + ylab("Likelihood") +
  geom_vline(xintercept=beta, col="blue")

fig2 <- arrangeGrob(alphaPP, betaPP, ncol=2, nrow=1, top="Baseline Model Profile Plot")
ggsave("../Figures/fig2.png", plot=fig2, width = 2 * 400, height = 1 * 400 + 20, units = "px", dpi=72)

# Figure 5 of the report
# Fetch data
trueAlphas <- numeric(nIterations*nSubjects)
trueBetas <- numeric(nIterations*nSubjects)
trueLLs <- numeric(nIterations*nSubjects)
fitAlphas <- numeric(nIterations*nSubjects)
fitBetas <- numeric(nIterations*nSubjects)
fitLLs <- numeric(nIterations*nSubjects)
k <- 1
for(i in 1:nIterations) {
  save_dir <- paste0("../Data/Simulated Data/Baseline Model Parameter Recovery/iteration-", i)
  for(j in 1:nSubjects) {
    # Load data
    load(paste0(save_dir, "/subject-", j, ".RData"))
    trueAlphas[k] <- trueAlpha
    trueBetas[k] <- trueBeta
    trueLLs[k] <- -sum(log(trueLL))
    fitAlphas[k] <- fitAlpha
    fitBetas[k] <- fitBeta
    fitLLs[k] <- -sum(log(fitLL))
    k <- k + 1
  }
}

# Plot
lmAlpha <- lm(fitAlphas ~ trueAlphas)
corAlpha <- cor(trueAlphas, fitAlphas, method="spearman")
alphaPlot <- ggplot(data.frame(trueAlphas=trueAlphas, fitAlphas=fitAlphas), aes(x=trueAlphas, y=fitAlphas)) +
  theme_classic() + ggtitle("(A) α recovery") +
  geom_point() + xlim(0, 1) + ylim(0, 1) + xlab("True α") + ylab("Estimated α") +
  geom_abline(aes(intercept=lmAlpha$coefficients[1], slope=lmAlpha$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=1, y=0, label=paste0("r = ", round(corAlpha, 2)))

lmBeta <- lm(fitBetas ~ trueBetas)
corBeta <- cor(trueBetas, fitBetas, method="spearman")
betaPlot <- ggplot(data.frame(trueBetas=trueBetas, fitBetas=fitBetas), aes(x=trueBetas, y=fitBetas)) +
  theme_classic() + ggtitle("(B) β recovery") +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab("True β") + ylab("Estimated β") +
  geom_abline(aes(intercept=lmBeta$coefficients[1], slope=lmBeta$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=5, y=0, label=paste0("r = ", round(corBeta, 2)))

lmLL <- lm(fitLLs ~ trueLLs)
corLL <- cor(trueLLs, fitLLs, method="spearman")
LLPlot <- ggplot(data.frame(trueLLs=trueLLs, fitLLs=fitLLs), aes(x=trueLLs, y=fitLLs)) +
  theme_classic() + ggtitle("(C) Likelihood recovery") +
  geom_point() + xlim(0, 75) + ylim(0, 75) + xlab("True likelihood") + ylab("Estimated likelihood") +
  geom_abline(aes(intercept=lmLL$coefficients[1], slope=lmLL$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=75, y=0, label=paste0("r = ", round(corLL, 2)))

fig5 <- arrangeGrob(alphaPlot, betaPlot, LLPlot, ncol=2, nrow=2, top="Baseline Model Parameter Recovery")
ggsave("../Figures/fig5.png", plot=fig5, width = 2 * 400, height = 2 * 400 + 20, units = "px", dpi=72)