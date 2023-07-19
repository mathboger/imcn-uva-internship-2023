# Script to use the simple delta rule + softmax (referred to as baseline model in the report) to analyze the experimental data
# Mathematical definitions and data wrangling

rm(list=ls())
library(rtdists)
library(DEoptim)
library(ggplot2)
library(gridExtra)

# Load the data
load("../Data/Original Data/data_exp2.RData")

# The variable 'dat' is the dataframe imported from data_exp2

### Global variables definitions

# Data structure
nParticipants <- 48
missingParticipant <- 38
unusefulParticipant <- 15 # All blocks have at least one NA
nBlocks <- 4
nTrials <- 128 # per block
nChoices <- 16 # 2 stimuli x 8 sets

# Setting parameters for optimization
DElower <- setNames(c(0.01, 0.01), 
                    c("alpha", "beta"))
DEupper <- setNames(c(1, 3),
                    c("alpha", "beta"))

# Convert chosen symbol into an integer to associate it to a VKF
choice_dictionary <- setNames(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
                              c("A", "B", "c", "e", "J", "K", "M", "n", "N", "O", "P", "s", "u", "W", "X", "y"))

### Auxiliary functions

# Helper function to map a vector of symbols to the corresponding integer according to the choice_dictionary
map_choices <- function(choice_symbol_list) {
  choices <- numeric(length(choice_symbol_list))
  for(i in 1:length(choice_symbol_list)) {
    choices[i] <- choice_dictionary[as.character(choice_symbol_list[i])]
  }
  return(choices)
}

# Function for fitting a model and calculating its likelihood (from vkf_lba.R)
ll_func <- function(pars, data, nTrials, nChoices, returnModel=FALSE) {
  
  # Fetch parameters from pars
  alpha <- pars["alpha"]
  beta <- pars["beta"]
  
  # Fetch outcome from data
  choice <- map_choices(data$choice_symbol)
  reward <- data$reward
  
  # Initialize variables
  Q <- matrix(NA, nrow=nTrials+1, ncol=nChoices)
  LL <- numeric(nTrials)
  
  Q[1,] <- rep(0, nChoices)
  
  for(i in 2:(nTrials+1)) {
    
    # Update variables (only change the Q-value associated to the choice made)
    Q[i,] <- Q[i-1,]
    Q[i, choice[i-1]] <- Q[i-1, choice[i-1]] + alpha*(reward[i] - Q[i-1, choice[i-1]])
    
    # Calculate likelihood
    pChoice <- exp(beta*Q[i-1,])/sum(exp(beta*Q[i-1,]))
    LL[i-1] <- pChoice[choice[i-1]]
  }
  
  if(returnModel) {
    model <- data.frame(Q=Q)
    return(model)
  }
  
  LL[LL == 0] <- 1e-20 # Adjustment for log in case we get a misfit
  return(LL)
}

# Wrapper function for DEoptim
wrapper <- function(pars, data, nTrials, nChoices) {
  pars <- setNames(pars, c("alpha", "beta"))
  LL <- ll_func(pars, data, nTrials, nChoices)
  return(-sum(log(LL)))
}

### Fit models (one per participant)

startTime <- Sys.time() # Keep track of how long it takes
for(participant in 4:nParticipants) {
  if(participant == missingParticipant) {
    next
  }
  if(participant == unusefulParticipant) {
    next
  }
  save_dir <- paste0("../Data/Fitted Models/Baseline Model")
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive=TRUE)
  reversalPoints <- numeric(nBlocks)
  for(block in 1:nBlocks) {
    data <- dat[dat$pp==participant & dat$block==block,]
    # In some participant x block combinations the data is incomplete (e.g., p=1, b=4)
    # Check if the data is complete or not before proceeding
    if(NA %in% data$rt) {
      reversalPoints[block] <- -1 # sentinel
      next
    }
    reversalPoints[block] <- data[data$trialNreversal == 0,]$TrialNumber
  }
  # Fit only using data from blocks without NAs
  data <- dat[dat$pp==participant & reversalPoints[dat$block] != -1,]
  out <- DEoptim(wrapper, DElower, DEupper, control=DEoptim.control(itermax=200) , data=data, nTrials=length(data$rt), nChoices=nChoices)
  pars <- setNames(out$optim$bestmem, c("alpha", "beta"))
  model <- ll_func(pars, data, length(data$rt), nChoices, returnModel=TRUE)
  ll <- ll_func(pars, data, length(data$rt), nChoices)
  # BIC; formally k*log(n) - 2*log(LL), where k is the number of parameters and n is the number of observations
  bic <- 8 * log(length(data$choice_symbol)) - 2 * log(-sum(log(ll)))
  save(reversalPoints, pars, model, ll, bic, data, file=paste0(save_dir, "/Participant-", participant, ".RData"))
}
endTime <- Sys.time()
print(endTime-startTime) # Check how long it took

### Data simulation from fit models

# Function to simulate experimental data from the fit model
generateData <- function(filename_participant_model) {
  
  # Fetch fit model
  load(filename_participant_model)
  nTrials <- length(data$rt)
  
  alpha <- pars["alpha"]
  beta <- pars["beta"]
  
  choice <- numeric(nTrials)
  reward <- numeric(nTrials)
  choseHigherReward <- numeric(nTrials) # needed to calculate accuracy
  
  # Initialize variables
  Q <- matrix(NA, nrow=nTrials+1, ncol=nChoices) # Q-value for RL
  
  Q[1,] <- rep(0, nChoices)
  
  # Adjustment for nTrials needed due to vector indexing
  for(i in 2:(nTrials+1)) {
    
    # Make observation
    # Get Q for the 2 stimuli shown in this trial and then reconvert to proper choice number
    temp_Q <- c(Q[i-1,choice_dictionary[data$stim_left[i-1]]], Q[i-1,choice_dictionary[data$stim_right[i-1]]])
    pChoice <- exp(beta*temp_Q)/sum(exp(beta*temp_Q))
    choice[i-1] <- rbinom(1, 1, pChoice[1])
    if(choice[i-1] == 1) {
      choice[i-1] <- choice_dictionary[data$stim_left[i-1]]
      pReward <- data$p_win_left[i-1] # reversal and block change already taken into consideration using this method
    } else {
      choice[i-1] <- choice_dictionary[data$stim_right[i-1]]
      pReward <- data$p_win_right[i-1]
    }
    reward[i-1] <- runif(1) < pReward
    if(pReward == data$p_win_correct[i-1]) {
      choseHigherReward[i-1] <- 1
    }
    
    # Update variables (only change the Q-values associated to the choice made)
    Q[i,] <- Q[i-1,]
    Q[i, choice[i-1]] <- Q[i-1, choice[i-1]] + alpha*(reward[i-1] - Q[i-1, choice[i-1]])
  }
  
  model <- data.frame(Q=Q)
  outcome <- data.frame(choice=choice, reward=reward, choseHigherReward=choseHigherReward)
  out <- setNames(list(model, outcome), c("model", "outcome"))
  return(out)
}

# Simulate data
# For RT, the average of the simulations is taken
# For accuracy, the mode of the simulations is taken
nSim <- 1000 # How many times the same trial should be simulated
sim_data <- list()
startTime <- Sys.time() # Keep track of how long it takes
for(participant in 1:nParticipants) {
  if(participant == missingParticipant) {
    next
  }
  if(participant == unusefulParticipant) {
    next
  }
  save_dir <- paste0("../Data/Simulated Data/Baseline Model Experiment Simulation")
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive=TRUE)
  for(i in 1:nSim) {
    print(paste0("Generating data for participant ", participant, ", sim ", i)) # Have some console output to see it's working
    sim_data[[i]] <- generateData(paste0("../Data/Fitted Models/Baseline Model/Participant-", participant, ".RData"))
    save(sim_data, file=paste0(save_dir, "/Participant-", participant, ".RData"))
  }
}
endTime <- Sys.time()
print(endTime-startTime) # Check how long it took

### Data retrieval and wrangling from files

# Set binning parameters for fetching data
nBins <- 10 # per block
bin_sizes <- numeric(nBins) + floor(nTrials/nBins)
# Distribute bin sizes as evenly as possible
i <- 1
while(sum(bin_sizes) < nTrials) {
  bin_sizes[i] <- bin_sizes[i] + 1
  i <- i + 1
  if(i > nBins) {
    i <- 1
  }
}

# Fetch data trial by trial from simulations
sim_agg_data <- list()
for(participant in 1:nParticipants) {
  if(participant == missingParticipant) {
    next
  }
  if(participant == unusefulParticipant) {
    next
  }
  load(paste0("../Data/Fitted Models/Baseline Model/Participant-", participant, ".RData"))
  load(paste0("../Data/Simulated Data/Baseline Model Experiment Simulation/Participant-", participant, ".RData"))
  overall_trial <- 1 # Go over simulated data trial by trial
  accuracy <- matrix(NA, nrow=nBlocks, ncol=nBins)
  for (block in 1:nBlocks) {
    # Check if the block was discarded
    if(reversalPoints[block] == -1) {
      next
    }
    # Fetch data bin by bin
    for(bin in 1:nBins) {
      choseHighestReward <- matrix(0, nrow=nSim, ncol=bin_sizes[bin])
      for(trial in 1:bin_sizes[bin]) {
        for(sim in 1:nSim) {
          choseHighestReward[sim, trial] <- sim_data[[sim]][["outcome"]]$choseHigherReward[overall_trial]
        }
        overall_trial <- overall_trial + 1
      }
      # Calculate accuracy for the bin
      sim_accuracy <- numeric(nSim)
      for(sim in 1:nSim) {
        sim_accuracy[sim] <- sum(choseHighestReward[sim,])/length(choseHighestReward[sim,])
      }
      accuracy[block, bin] <- mean(sim_accuracy)
    }
  }
  participant_data <- setNames(list(accuracy), c("accuracy"))
  sim_agg_data[[participant]] <- participant_data
}

# Collapse fetched data through participants
nParticipants_per_block <- numeric(nBlocks)
sim_accuracy <- matrix(NA, nrow=nBlocks, ncol=nBins) # Final variable for plot
for(block in 1:nBlocks) {
  for(bin in 1:nBins) {
    data_accuracy <- rep(-1, nParticipants)
    for(participant in 1:nParticipants) {
      if(participant == missingParticipant) {
        next
      }
      if(participant == unusefulParticipant) {
        next
      }
      if(is.null(sim_agg_data[[participant]])) {
        next
      }
      if(!is.na(sim_agg_data[[participant]]$accuracy[block, bin])) {
        data_accuracy[participant] <- sim_agg_data[[participant]]$accuracy[block, bin]
      }
    }
    nParticipants_per_block[block] <- length(data_accuracy[data_accuracy != -1])
    sim_accuracy[block, bin] <- mean(data_accuracy[data_accuracy != -1]) # Discard unavailable data
  }
}

# Fetch real data
real_agg_data <- list()
for(participant in 1:nParticipants) {
  if(participant == missingParticipant) {
    next
  }
  if(participant == unusefulParticipant) {
    next
  }
  load(paste0("../Data/Fitted Models/Baseline Model/Participant-", participant, ".RData")) # For checking reversal point sentinel
  accuracy <- matrix(NA, nrow=nBlocks, ncol=nBins)
  for (block in 1:nBlocks) {
    overall_trial <- 1 # Go over simulated data trial by trial; must be reseted each block due to formatting of the "dat" data frame
    # Check if the block was discarded
    if(reversalPoints[block] == -1) {
      next
    }
    # Fetch data bin by bin
    for(bin in 1:nBins) {
      choseHighestReward <- numeric(bin_sizes[bin])
      for(trial in 1:bin_sizes[bin]) {
        choseHighestReward[trial] <- dat[dat$pp == participant & dat$block == block & dat$TrialNumber == overall_trial,]$choiceIsHighP
        overall_trial <- overall_trial + 1
      }
      # Calculate accuracy for the bin
      real_accuracy <- sum(choseHighestReward)/length(choseHighestReward)
      accuracy[block, bin] <- real_accuracy
    }
  }
  participant_data <- setNames(list(accuracy), c("accuracy"))
  real_agg_data[[participant]] <- participant_data
}

# Collapse fetched data through participants
real_accuracy <- matrix(NA, nrow=nBlocks, ncol=nBins) # Final variable for plot
real_accuracy_std_error <- matrix(NA, nrow=nBlocks, ncol=nBins) # Final variable for plot
for(block in 1:nBlocks) {
  for(bin in 1:nBins) {
    data_accuracy <- rep(-1, nParticipants)
    for(participant in 1:nParticipants) {
      if(participant == missingParticipant) {
        next
      }
      if(participant == unusefulParticipant) {
        next
      }
      if(is.null(real_agg_data[[participant]])) {
        next
      }
      if(!is.na(real_agg_data[[participant]]$accuracy[block, bin])) {
        data_accuracy[participant] <- real_agg_data[[participant]]$accuracy[block, bin]
      }
    }
    real_accuracy[block, bin] <- mean(data_accuracy[data_accuracy != -1]) # Discard unavailable data
    real_accuracy_std_error[block, bin] <- sd(data_accuracy[data_accuracy != -1]) / sqrt(nParticipants_per_block[block])
  }
}

### Plotting

# Figure 8 of the report
bins <- seq(1, nBins, 1) # for x axis
plots <- list()
for(block in 1:nBlocks) {
  
  # Accuracy
  plot <- ggplot(data.frame(sim_accuracy=sim_accuracy[block,], real_accuracy=real_accuracy[block,], real_accuracy_std_error=real_accuracy_std_error[block,], bins=bins)) +
    theme_classic() + ylim(0, 1) + ggtitle(paste0("Block ", block)) +
    scale_x_continuous(breaks=bins, labels=bins) + xlab("Bin") + ylab("Accuracy") +
    geom_line(aes(x=bins, y=sim_accuracy), color="blue") +
    geom_point(aes(x=bins, y=sim_accuracy), color="blue") +
    geom_point(aes(x=bins, y=real_accuracy)) +
    geom_line(aes(x=bins, y=real_accuracy)) +
    geom_ribbon(aes(x=bins, ymin=(real_accuracy - real_accuracy_std_error), ymax=(real_accuracy + real_accuracy_std_error)), alpha=0.2)
  
  plots[[block]] <- plot
}

fig8 <- arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                     ncol=4, nrow=1, top="Baseline Model Fit to Data")
ggsave("../Figures/fig8.png", plot=fig8, width = 4 * 400, height = 1 * 400 + 20, units = "px", dpi=72)