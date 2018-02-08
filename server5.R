# MCMC estimating stochastic processes
# 3 unknowns: beta, gamma, and the time of infection for each individual
# 17/11/17

###############
## Read data ##
###############
# setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")
run_stoch <- read.csv("data_pop100_b2e-2_g7e-2_15.csv")

###########
## Input ##
###########

# Guesses for beta and gamma
beta = 0.005
gamma = 0.005

# Proposal function SDs
prop_sd_beta = 0.002
prop_sd_gamma = 0.007

# inf_period <- 10 # Guess an infectious period for infectious curve starting point
inf_period <- ceiling(1/gamma) # mean infectious period calculated from gamma

# Number of runs
iterations = 3500000 # How many iterations MCMC is running for
divisor = 1000 # How often runs are being saved

#############################
## Approximate new_I and I ##
#############################

N = run_stoch$S[1] + run_stoch$I[1] + run_stoch$R[1]

timestep <- run_stoch$time[2] - run_stoch$time[1]

inf_timestep <- inf_period/timestep # translates infectious days into no. of timesteps
epi_time <- max(run_stoch$time) + inf_period
times <- seq(0, epi_time, by = timestep)

zero_array <- array(0, dim = c(inf_timestep, ncol(run_stoch)))
colnames(zero_array) <- c("time","S", "I", "R", "new_I", "new_R")
run_stoch <- rbind(zero_array, run_stoch)
run_stoch$time <- times

for (i in 1:nrow(run_stoch)){
  run_stoch$guess_new_I[i] <- {
    run_stoch$new_R[i+(inf_timestep)] # guess newly infected, translates days into no. of timesteps
  }
  
  run_stoch$guess_I[i] <- if(i == 1){
    run_stoch$guess_new_I[i] 
  } else{
    run_stoch$guess_I[i-1] + run_stoch$guess_new_I[i] - run_stoch$new_R[i]
  }
  run_stoch$guess_S[i] <- N - (run_stoch$R[i] + run_stoch$guess_I[i])
}

run_stoch[is.na(run_stoch)] <- 0

##################################
## Likelihood, prior, posterior ##
##################################

# Likelihood distribution
bg_likelihood <- function(param){
  beta = as.numeric(param[1,1,1])
  gamma = as.numeric(param[2,1,1])
  I = as.numeric(param[,1,2]) # Infected
  new_I = as.numeric(param[,1,3]) # Newly infected
  
  total = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:nrow(run_stoch -1)){
    S = (N - I[i] - run_stoch$R[i]) # Susceptibles
    betalikelihood = dbinom(new_I[i+1], S, (1-(exp(-beta*I[i]*timestep))), log = T)
    gammalikelihood = dbinom(run_stoch$new_R[i+1], I[i], (1-(exp(-gamma*timestep))), log = T)
    
    # To deal with -Inf (i.e. log(0))
    if (is.na(new_I[i+1]) == F & betalikelihood == -Inf){
      betalikelihood = -1000
    }
    if (is.na(new_I[i+1]) == F & gammalikelihood == -Inf){
      gammalikelihood = -1000
    }
    total[i] = (betalikelihood + gammalikelihood)
  }
  
  ll = sum(total, na.rm = T)
  return(ll)
}

inf_likelihood <- function(param){
  beta = as.numeric(param[1,2,1])
  gamma = as.numeric(param[2,2,1])
  I = as.numeric(param[,1,2]) # Infected
  new_I = as.numeric(param[,1,3]) # Newly infected
  
  total = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:(nrow(run_stoch) -1)){
    S = (N - I[i] - run_stoch$R[i]) # Susceptibles
    betalikelihood = dbinom(new_I[i+1], S, (1-(exp(-beta*I[i]*timestep))), log = T)
    gammalikelihood = dbinom(run_stoch$new_R[i+1], I[i], (1-(exp(-gamma*timestep))), log = T)
    
    # To deal with -Inf (i.e. log(0))
    if (is.na(new_I[i+1]) == F & betalikelihood == -Inf){
      betalikelihood = -1000
    }
    if (is.na(new_I[i+1]) == F & gammalikelihood == -Inf){
      gammalikelihood = -1000
    }
    total[i] = (betalikelihood + gammalikelihood)
  }
  
  ll = sum(total, na.rm = T)
  return(ll)
}

# Prior distribution
prior <- function(param){
  beta = as.numeric(param[1,1,1])
  gamma = as.numeric(param[2,1,1])
  
  betaprior = dunif(beta, min = 0, max = 100, log = T)
  gammaprior = dunif(gamma, min = 0, max = 100, log = T)
  
  return(betaprior + gammaprior)
}

# Posterior distribution
inf_posterior <- function(param){
  return (inf_likelihood(param) + prior(param))
}

bg_posterior <- function(param){
  return (bg_likelihood(param) + prior(param))
}

# Proposal function for infectious
inf_proposalfunction <- function(param){
  changed_I <- sample(nrow(run_stoch), 1)
  inf_list <- c(-1, 1) # used for choosing -1 or +1 randomly
  
  inf <- sample(c(-1, 1), 1) # will the change at that timepoint be + or - n I
  
  neighbour <- if (changed_I == 1){
    changed_I + 1
  } else if (changed_I == nrow(run_stoch)){
    changed_I - 1
  } else{
    changed_I + sample(inf_list, 1) # will choose which neighbouring timepoint is also affected
  }
  
  param[changed_I, 1, 3] = param[changed_I, 1, 3] + inf
  param[neighbour, 1, 3] = param[neighbour, 1, 3] - inf
  
  # Calculate total infectious curve from changed new_I values
  for (i in 1:nrow(run_stoch)){ 
    param[i,1,2] <- if(i == 1){
      param[i,1,3] 
    } else{
      param[i-1,1,2] + param[i,1,3] - run_stoch$new_R[i]
    }
  }
  return(param)
}

# Proposal function for beta and gamma
proposalfunction <- function(param){
  param[1, 1, 1] = rnorm(1, mean = param[1, 1, 1], sd = prop_sd_beta) # beta proposal
  param[2, 1, 1] = rnorm(1, mean = param[2, 1, 1], sd = prop_sd_gamma) # gamma proposal
  return(param)
}

##########
## MCMC ##
##########

run_metropolis_MCMC <- function(startvalue, iterations){
  divisor = divisor # the interval at which chain values are saved
  chain = array(dim = c(nrow(startvalue), (iterations/divisor), ncol(startvalue))) # Array for storing chain data
  temp_chain = array(dim = c(nrow(startvalue), 2, ncol(startvalue))) # Temporary array used for single iterations
  
  chain[1,1,1] = startvalue[1,1] # beta
  chain[2,1,1] = startvalue[2,1] # gamma
  chain[,1,2] = startvalue[,2] # all infectious
  chain[,1,3] = startvalue[,3] # newly infectious
  
  temp_chain[,1,] = chain[,1,]
  
  for (i in 1:iterations){
    
    # Beta and gamma
    proposal = proposalfunction(temp_chain[,,]) # beta proposal and gamma proposal
    
    if (proposal[1,1,1] < 0.0 & proposal[2,1,1] < 0.0 ){
      temp_chain[1:2,2,1] = temp_chain[1:2,1,1]
      
    } else if (proposal[1,1,1] < 0.0 & proposal[2,1,1] >= 0.0){
      proposal[1,1,1] = temp_chain[1,1,1]
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[,,])
      if (log(runif(1)) < probab){
        temp_chain[,2,1] = proposal[,1,1]
      } else{
        temp_chain[,2,1] = temp_chain[,1,1]
      }
      
    } else if (proposal[1,1,1] >= 0.0 & proposal[2,1,1] < 0.0){
      proposal[2,1,1] = temp_chain[2,1,1]
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[,,])
      if (log(runif(1)) < probab){
        temp_chain[,2,1] = proposal[,1,1]
      } else{
        temp_chain[,2,1] = temp_chain[,1,1]
      }
      
    } else{
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[,,])
      if (log(runif(1)) < probab){
        temp_chain[,2,1] = proposal[,1,1]
      } else{
        temp_chain[,2,1] = temp_chain[,1,1]
      }
    }
    
    # Infectious
    inf_proposal = inf_proposalfunction(temp_chain[,,])
    
    # Check susceptibles for negative numbers
    S = array(dim = c(nrow(startvalue)))
    for (time in (1:nrow(startvalue))){
      S[time] = (N - (inf_proposal[time,1,2] + run_stoch$R[time]))
    }
    
    if(min(inf_proposal[,1,3]) < 0 | min(S) < 0){
      temp_chain[,2,2:3] = temp_chain[,1,2:3]
    } else{
      inf_probab = inf_posterior(inf_proposal) - inf_posterior(temp_chain[,,])
      if (log(runif(1)) < inf_probab){
        temp_chain[,2,2:3] = inf_proposal[,1,2:3]
      } else{
        temp_chain[,2,2:3] = temp_chain[,1,2:3]
      }
    }
    
    # Save the iteration if it is divisible by the divisor without residuals
    if (i%%divisor == 0) {
      chain[,(i/divisor),] = temp_chain[,2,]
      chain[3,(i/divisor),1] = inf_likelihood(temp_chain[,,]) # saving the log-likelihood
    }    
    
    # Re-set temporary chain for next iteration
    temp_chain[,1,] = temp_chain[,2,]
    
    # Print every nth iteration, to know how far along run is
    if (i%%(iterations/20) == 0 | i == 1) {
      print(i)
    }
  }
  return(chain)
}

# Where to start the chain for beta, gamma, and I
startvalue <- array(dim = c(nrow(run_stoch), 3))
startvalue[1,1] <- beta # beta guess
startvalue[2,1] <- gamma # gamma guess
startvalue[,2] <- run_stoch$guess_I # I guess 
startvalue[,3] <- run_stoch$guess_new_I # new I guess

# Run the MCMC
set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)

########################################################################################################################

########################################################################################################################

#################
## Saving data ##
#################

# Prepare data for saving
## Saving beta, gamma, and log likelihood
beta_gamma_loglik <- t(chain[1:3,,1]) # transpose beta, gamma, and log likelihood data so that rows become columns etc.
colnames(beta_gamma_loglik) <- c("beta", "gamma", "loglik") # name the columns
## Saving infectious data
inf_data <- chain[,,2]
timeframe <- array(dim = nrow(inf_data))

for (i in 1:nrow(inf_data)){
  timeframe[i] = i * timestep 
}
inf_data2 <- cbind(timeframe, inf_data)

# setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")

# Beta, gamma, and likelihood data
write.csv(data.frame(beta_gamma_loglik), file = "mcmc_pop100_b2e-2_g7e-2_15_loglik.csv", row.names = FALSE)

# Infectious curve data
write.csv(data.frame(inf_data2), file = "mcmc_pop100_b2e-2_g7e-2_15_infectious.csv", row.names = FALSE)
