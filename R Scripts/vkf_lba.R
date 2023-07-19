# Script to test the VKF + LBA model (referred as target model in the report)
# Mathematical definitions and parameter recovery (data & plots)

rm(list=ls())
library(rtdists)
library(DEoptim)
library(ggplot2)
library(gridExtra)

### Generate data following the predictive binary model of the binary VKF
### Connect to LBA with M (M*V) (i.e., M serves as a Q-value)

# VKF parameters
m0 <- 0.5 # Initial state
v0 <- 1 # Initial volatility
lambda <- 0.5 # Volatility update rate
omega <- 1 # Observation noise

# LBA parameters
t0 <- .2 # Non-decision time
dr <- c(1, 1) # Mean drift rates for stimulus 1 and 2 (avoided v not to confound with volatility)
B <- 1 # Threshold
A <- .5 # Between trial variability in start point
s <- c(1,1) # Between trial variability in drift rate; fixed
b <- B + A

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

# Function that takes in VKF, LBA, and experiment parameters and outputs simulated data along the model
generateData <- function(nTrials, nChoices, reversalPoint, pReward, m0, v0, lambda, omega, t0, dr, A, B) {
  
  # Finalize LBA parameters
  b <- B + A
  s <- c(1, 1)
  
  # Initialize variables
  # VKF
  K <- matrix(NA, nrow=nTrials, ncol=nChoices)
  Alpha <- matrix(NA, nrow=nTrials, ncol=nChoices)
  M <- matrix(NA, nrow=nTrials, ncol=nChoices)
  W <- matrix(NA, nrow=nTrials, ncol=nChoices)
  CoW <- matrix(NA, nrow=nTrials, ncol=nChoices)
  V <- matrix(NA, nrow=nTrials, ncol=nChoices)
  # Output
  choice <- numeric(nTrials)
  RT <- numeric(nTrials)
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
    trial <- rLBA(1, A=A, b=b, t0=t0, mean_v=dr*M[i-1,], sd_v=s) # LBA
    choice[i] <- trial$response
    RT[i] <- trial$rt
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
  
  model <- data.frame(K=K, Alpha=Alpha, M=M, W=W, CoW=CoW, V=V, choice=choice, RT=RT, reward=reward)
  return(model)
}

# Usage example
model <- generateData(nTrials, nChoices, reversalPoint, pReward, m0, v0, lambda, omega, t0, dr, A, B)

### Likelihood

# Function that takes in the observed outcome and the proposed model parameters and outputs 
# the likelihood of that model having generated the data or the model itself (intermediate variables)
ll_func <- function(RT, choice, reward, nTrials, nChoices, m0, v0, lambda, omega, t0, dr, A, b, returnModel=FALSE) {
  
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
  LL <- dLBA(RT[2:nTrials], choice[2:nTrials], A, b, t0, mean_v=list(dr[1]*M[,1], dr[2]*M[,2]), sd_v=list(numeric(nTrials) + 1, numeric(nTrials) + 1), silent=TRUE)
  LL <- c(1, LL) # Placeholder that won't change the value of sum(log()) but will make LL the right size to couple it to a model
  LL[LL == 0] <- 1e-20 # Adjustment for log in case we get a misfit
  
  if(returnModel) {
    model <- data.frame(K=K, Alpha=Alpha, M=M, W=W, CoW=CoW, V=V, choice=choice, RT=RT, reward=reward)
    return(model)
  }
  return(LL)
}

### Profile Data Generation

rangeM0 <- seq(0.01, 5, 0.01)
rangeV0 <- seq(0.01, 5, 0.01)
rangeLambda <- seq(0.01, 0.99, 0.01)
rangeOmega <- seq(0.01, 5, 0.01)
rangeT0 <- seq(0.2, 1, 0.01)
rangeDR <- seq(0.01, 5, 0.01)
rangeA <- seq(0.01, 5, 0.01)
rangeB <- seq(0.01, 5, 0.01)

LLM0 <- numeric(length(rangeM0))
for(i in 1:length(rangeM0)) {
  LLM0[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, rangeM0[i], v0, lambda, omega, t0, dr, A, A+B)))
}
LLV0 <- numeric(length(rangeV0))
for(i in 1:length(rangeV0)) {
  LLV0[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, m0, rangeV0[i], lambda, omega, t0, dr, A, A+B)))
}
LLLambda <- numeric(length(rangeLambda))
for(i in 1:length(rangeLambda)) {
  LLLambda[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, m0, v0, rangeLambda[i], omega, t0, dr, A, A+B)))
}
LLOmega <- numeric(length(rangeOmega))
for(i in 1:length(rangeOmega)) {
  LLOmega[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, m0, v0, lambda, rangeOmega[i], t0, dr, A, A+B)))
}
LLT0 <- numeric(length(rangeT0))
for(i in 1:length(rangeT0)) {
  LLT0[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, m0, v0, lambda, omega, rangeT0[i], dr, A, A+B)))
}
LLDR <- numeric(length(rangeDR))
for(i in 1:length(rangeDR)) {
  LLDR[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, m0, v0, lambda, omega, t0, c(rangeDR[i], rangeDR[i]), A, A+B)))
}
LLA <- numeric(length(rangeA))
for(i in 1:length(rangeA)) {
  LLA[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, m0, v0, lambda, omega, t0, dr, rangeA[i], rangeA[i]+B)))
}
LLB <- numeric(length(rangeB))
for(i in 1:length(rangeB)) {
  LLB[i] <- sum(log(ll_func(model$RT, model$choice, model$reward, nTrials, nChoices, m0, v0, lambda, omega, t0, dr, A, A+rangeB[i])))
}
save(rangeM0, rangeV0, rangeLambda, rangeOmega, rangeT0, rangeDR, rangeA, rangeB, LLM0, LLV0, LLLambda, LLOmega, LLT0, LLDR, LLA, LLB, m0, v0, lambda, omega, t0, dr, A, B, file="../Data/Simulated Data/Target Model Profile/profile.RData")

### Parameter Recovery Data Generation

# Wrapper function for DEoptim
wrapper <- function(pars, nTrials, nChoices, trueModel) {
  LL <- ll_func(trueModel$RT, trueModel$choice, trueModel$reward, nTrials, nChoices, pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7], pars[7] + pars[8])
  return(-sum(log(LL)))
}

startTime <- Sys.time() # Keep track of how long it takes
for(i in 1:nIterations) {
  # Get save file / directory
  save_dir = paste0("../Data/Simulated Data/Target Model Parameter Recovery/iteration-", i)
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive=TRUE)
  
  # Sample new ground truth
  trueM0 <- runif(1, 0.01, 5)
  trueV0 <- runif(1, 0.01, 5)
  trueLambda <- runif(1, 0.01, 0.99)
  trueOmega <- runif(1, 0.01, 5)
  trueT0 <- runif(1, 0.2, 1.5)
  trueDR <- runif(1, 0.01, 5)
  trueA <- runif(1, 0.01, 5)
  trueB <- runif(1, 0.01, 5)
  
  for(j in 1:nSubjects) {
    
    # Get true model and likelihood
    trueModel <- generateData(nTrials, nChoices, reversalPoint, pReward, trueM0, trueV0, trueLambda, trueOmega, trueT0, trueDR, trueA, trueB)
    trueLL <- ll_func(trueModel$RT, trueModel$choice, trueModel$reward, nTrials, nChoices, trueM0, trueV0, trueLambda, trueOmega, trueT0, trueDR, trueA, trueA + trueB)
    
    # Get fit model and likelihood
    out <- DEoptim(wrapper, lower=c(0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.01, 0.01), upper=c(5, 5, 0.99, 5, 1.5, 5, 5, 5), nTrials=nTrials, nChoices=nChoices, trueModel=trueModel)
    fitM0 <- out$optim$bestmem[1]
    fitV0 <- out$optim$bestmem[2]
    fitLambda <- out$optim$bestmem[3]
    fitOmega <- out$optim$bestmem[4]
    fitT0 <- out$optim$bestmem[5]
    fitDR <- out$optim$bestmem[6]
    fitA <- out$optim$bestmem[7]
    fitB <- out$optim$bestmem[8]
    fitModel <- ll_func(trueModel$RT, trueModel$choice, trueModel$reward, nTrials, nChoices, fitM0, fitV0, fitLambda, fitOmega, fitT0, fitDR, fitA, fitA + fitB, returnModel=TRUE)
    fitLL <- ll_func(trueModel$RT, trueModel$choice, trueModel$reward, nTrials, nChoices, fitM0, fitV0, fitLambda, fitOmega, fitT0, fitDR, fitA, fitA + fitB)
    
    # Save ground truth and fit to save file
    save(trueM0, trueV0, trueLambda, trueOmega, trueT0, trueDR, trueA, trueB, 
         fitM0, fitV0, fitLambda, fitOmega, fitT0, fitDR, fitA, fitB, 
         trueModel, fitModel, trueLL, fitLL, file=paste0(save_dir, "/subject-", j, ".RData"))
  }
}
endTime <- Sys.time()
print(endTime-startTime) # Check how long it took

### Plotting

# Figure 4 of the report
# Fetch data
load("../Data/Simulated Data/Target Model Profile/profile.RData")

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

APP <- ggplot(data.frame(rangeA=rangeA, LLA=LLA), aes(x=rangeA, y=LLA)) +
  theme_classic() + ggtitle("(E) A") +
  geom_line() + xlab("A") + ylab("Likelihood") +
  geom_vline(xintercept=A, col="blue")

BPP <- ggplot(data.frame(rangeB=rangeB, LLB=LLB), aes(x=rangeB, y=LLB)) +
  theme_classic() + ggtitle("(F) B") +
  geom_line() + xlab("B") + ylab("Likelihood") +
  geom_vline(xintercept=B, col="blue")

DRPP <- ggplot(data.frame(rangeDR=rangeDR, LLDR=LLDR), aes(x=rangeDR, y=LLDR)) +
  theme_classic() + ggtitle("(G) υ") +
  geom_line() + xlab("υ") + ylab("Likelihood") +
  geom_vline(xintercept=dr[1], col="blue")

t0PP <- ggplot(data.frame(rangeT0=rangeT0, LLT0=LLT0), aes(x=rangeT0, y=LLT0)) +
  theme_classic() + ggtitle(expression((H)~t[0])) +
  geom_line() + xlab(expression(t[0])) + ylab("Likelihood") +
  geom_vline(xintercept=t0, col="blue")

fig4 <- arrangeGrob(m0PP, v0PP, lambdaPP, omegaPP, APP, BPP, DRPP, t0PP, ncol=2, nrow=4, top="Target Model Profile Plot")
ggsave("../Figures/fig4.png", plot=fig4, width = 2 * 400, height = 4 * 400 + 20, units = "px", dpi=72)

# Figure 7 of the report
# Fetch data
trueM0s <- numeric(nIterations*nSubjects)
trueV0s <- numeric(nIterations*nSubjects)
trueLambdas <- numeric(nIterations*nSubjects)
trueOmegas <- numeric(nIterations*nSubjects)
trueT0s <- numeric(nIterations*nSubjects)
trueDRs <- numeric(nIterations*nSubjects)
trueAs <- numeric(nIterations*nSubjects)
trueBs <- numeric(nIterations*nSubjects)
trueLLs <- numeric(nIterations*nSubjects)
fitM0s <- numeric(nIterations*nSubjects)
fitV0s <- numeric(nIterations*nSubjects)
fitLambdas <- numeric(nIterations*nSubjects)
fitOmegas <- numeric(nIterations*nSubjects)
fitT0s <- numeric(nIterations*nSubjects)
fitDRs <- numeric(nIterations*nSubjects)
fitAs <- numeric(nIterations*nSubjects)
fitBs <- numeric(nIterations*nSubjects)
fitLLs <- numeric(nIterations*nSubjects)
k <- 1
for(i in 1:nIterations) {
  save_dir <- paste0("../Data/Simulated Data/Target Model Parameter Recovery/iteration-", i)
  for(j in 1:nSubjects) {
    # Load data
    load(paste0(save_dir, "/subject-", j, ".RData"))
    trueM0s[k] <- trueM0
    trueV0s[k] <- trueV0
    trueLambdas[k] <- trueLambda
    trueOmegas[k] <- trueOmega
    trueT0s[k] <- trueT0
    trueDRs[k] <- trueDR
    trueAs[k] <- trueA
    trueBs[k] <- trueB
    trueLLs[k] <- -sum(log(trueLL))
    fitM0s[k] <- fitM0
    fitV0s[k] <- fitV0
    fitLambdas[k] <- fitLambda
    fitOmegas[k] <- fitOmega
    fitT0s[k] <- fitT0
    fitDRs[k] <- fitDR
    fitAs[k] <- fitA
    fitBs[k] <- fitB
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

lmA <- lm(fitAs ~ trueAs)
corA <- cor(trueAs, fitAs, method="spearman")
APlot <- ggplot(data.frame(trueAs=trueAs, fitAs=fitAs), aes(x=trueAs, y=fitAs)) +
  theme_classic() + ggtitle("(E) A recovery") +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab("True A") + ylab("Estimated A") +
  geom_abline(aes(intercept=lmA$coefficients[1], slope=lmA$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=5, y=0, label=paste0("r = ", round(corA, 2)))

lmB <- lm(fitBs ~ trueBs)
corB <- cor(trueBs, fitBs, method="spearman")
BPlot <- ggplot(data.frame(trueBs=trueBs, fitBs=fitBs), aes(x=trueBs, y=fitBs)) +
  theme_classic() + ggtitle("(F) B recovery") +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab("True B") + ylab("Estimated B") +
  geom_abline(aes(intercept=lmB$coefficients[1], slope=lmB$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=5, y=0, label=paste0("r = ", round(corB, 2)))

lmDR <- lm(fitDRs ~ trueDRs)
corDR <- cor(trueDRs, fitDRs, method="spearman")
DRPlot <- ggplot(data.frame(trueDRs=trueDRs, fitDRs=fitDRs), aes(x=trueDRs, y=fitDRs)) +
  theme_classic() + ggtitle("(G) υ recovery") +
  geom_point() + xlim(0, 5) + ylim(0, 5) + xlab("True υ") + ylab("Estimated υ") +
  geom_abline(aes(intercept=lmDR$coefficients[1], slope=lmDR$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=5, y=0, label=paste0("r = ", round(corDR, 2)))

lmT0 <- lm(fitT0s ~ trueT0s)
corT0 <- cor(trueT0s, fitT0s, method="spearman")
t0Plot <- ggplot(data.frame(trueT0s=trueT0s, fitT0s=fitT0s), aes(x=trueT0s, y=fitT0s)) +
  theme_classic() + ggtitle(expression((H)~t[0]~recovery)) +
  geom_point() + xlim(0, 1.5) + ylim(0, 1.5) + xlab(expression(True~t[0])) + ylab(expression(Estimated~t[0])) +
  geom_abline(aes(intercept=lmT0$coefficients[1], slope=lmT0$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=1.5, y=0, label=paste0("r = ", round(corT0, 2)))

lmLL <- lm(fitLLs ~ trueLLs)
corLL <- cor(trueLLs, fitLLs, method="spearman")
LLPlot <- ggplot(data.frame(trueLLs=trueLLs, fitLLs=fitLLs), aes(x=trueLLs, y=fitLLs)) +
  theme_classic() + ggtitle("(I) Likelihood recovery") +
  geom_point() + xlim(1000, 2750) + ylim(1000, 2750) + xlab("True likelihood") + ylab("Estimated likelihood") +
  geom_abline(aes(intercept=lmLL$coefficients[1], slope=lmLL$coefficients[2])) +
  geom_abline(aes(intercept=0, slope=1), color="blue") +
  geom_text(x=2750, y=1000, label=paste0("r = ", round(corLL, 2)))

fig7 <- arrangeGrob(m0Plot, v0Plot, lambdaPlot, omegaPlot, APlot, BPlot, DRPlot, t0Plot, LLPlot, ncol=2, nrow=5, top="Target Model Parameter Recovery")
ggsave("../Figures/fig7.png", plot=fig7, width = 2 * 400, height = 5 * 400 + 20, units = "px", dpi=72)