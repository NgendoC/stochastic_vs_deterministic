# MCMC estimating stochastic processes with a deterministic likelihood
# 3 unknowns: beta, gamma, and the time of infection for each individual
# Assuming that time of infection is deterministic
# 18/12/17

#######################
## Required packages ##
#######################
#install.packages("deSolve")
#if (!require("deSolve")) install.packages("deSolve")
library("deSolve") #package for solving differential equations

###############
## Read data ##
###############
# setwd("/home/evelina/Development/stochastic_vs_deterministic")
run_stoch <- read.csv("data_pop50_R6_g0.15_2.csv")

###########
## Input ##
###########

# Initial beta and gamma guesses
N = run_stoch$S[1] + run_stoch$I[1] + run_stoch$R[1]

R0 <-7.2 # R0 = beta*N/gamma
gamma <- 0.1
beta <- R0*gamma/N

# Number of runs
iterations =  2500000 # how many iterations the MCMC will be running for
divisor = 1000 # how often runs are being saved

# Proposal function SDs
prop_sd_beta = beta/4
prop_sd_gamma = gamma/4

##################################
## Likelihood, prior, posterior ##
##################################

# Time
timestep <- run_stoch$time[2] - run_stoch$time[1]
end <- max(run_stoch$time)
times <- seq(0, end, by = timestep)

# Initial population: N-1 susceptible, 1 infectious, 0 recovered
init.values = c(
  S = run_stoch$S[1],
  I = run_stoch$I[1],
  R = run_stoch$R[1]
)
N = sum(init.values)

sir <- function(time, state, param) {
  
  # define model parameters in term of the natural parameters
  beta <- param[1,1] 
  gamma <- param[2,1]
  
  with(as.list(c(state, param)), {
    
    dS <- -beta * S * I 
    dI <- (beta * S * I) -(gamma * I)
    dR <-  gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Likelihood distribution
bg_likelihood <- function(param){
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  
  det_sir <- ode(y = init.values, times = times, func = sir, parms = param)
  det_sir <- as.data.frame(det_sir)
  
  total = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:nrow(run_stoch -1)){
    I = round(det_sir$I[i])
    S = (N - (I + run_stoch$R[i])) # Susceptibles for timestep i
    new_I = if (i == 1){
      I
    } else {
      (round(det_sir$I[i+1]) - I + run_stoch$R[i+1] - run_stoch$R[i]) # new I for timestep i+1
    }
    
    if (is.na(new_I) == F & (S<0 | new_I<0 | I<0)){
      betalikelihood = -1000
      gammalikelihood = -1000
    } else {
      betalikelihood = dbinom(new_I, S, (1-(exp(-beta*I*timestep))), log = T)
      gammalikelihood = dbinom(run_stoch$new_R[i+1], I, (1-(exp(-gamma*timestep))), log = T)
    }
    
    # To deal with -Inf (i.e. log(0))
    if (is.na(new_I) == F & betalikelihood == -Inf){
      betalikelihood = -1000
    }
    if (is.na(new_I) == F & gammalikelihood == -Inf){
      gammalikelihood = -1000
    }
    total[i] = (betalikelihood + gammalikelihood)
  }
  
  ll = sum(total, na.rm = T)
  return(ll)
}

# Prior distribution
prior <- function(param){
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  
  betaprior = dunif(beta, min = 0, max = 100, log = T)
  gammaprior = dunif(gamma, min = 0, max = 100, log = T)
  
  return(betaprior + gammaprior)
}

# Posterior distribution
bg_posterior <- function(param){
  return (bg_likelihood(param) + prior(param))
}

# Proposal function for beta and gamma
proposalfunction <- function(param){
  param[1, 1] = rnorm(1, mean = param[1, 1], sd = prop_sd_beta) # beta proposal
  param[2, 1] = rnorm(1, mean = param[2, 1], sd = prop_sd_gamma) # gamma proposal
  return(param)
}

##########
## MCMC ##
##########

run_metropolis_MCMC <- function(startvalue, iterations){
  divisor = divisor # the interval at which chain values are saved
  chain = array(dim = c((nrow(startvalue)+1), (iterations/divisor))) # Array for storing chain data
  temp_chain = array(dim = c((nrow(startvalue)+1), 2)) # Temporary array used for single iterations
  
  chain[1,1] = startvalue[1] # beta
  chain[2,1] = startvalue[2] # gamma
  chain[3,1] = 0             # likelihood
  
  
  temp_chain[,1] = chain[,1]
  
  for (i in 1:iterations){
    
    # Beta and gamma
    proposal = proposalfunction(temp_chain[1:2,]) # beta proposal and gamma proposal
    
    if (proposal[1,1] < 0.0 & proposal[2,1] < 0.0 ){
      temp_chain[1:2,2] = temp_chain[1:2,1]
      temp_chain[3,2] = temp_chain[3,1]
      
    } else if (proposal[1,1] < 0.0 & proposal[2,1] >= 0.0){
      proposal[1,1] = temp_chain[1,1]
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[1:2,])
      if (log(runif(1)) < probab){
        temp_chain[1:2,2] = proposal[1:2,1]
        temp_chain[3,2] = bg_likelihood(proposal)
      } else{
        temp_chain[1:2,2] = temp_chain[1:2,1]
        temp_chain[3,2] = bg_likelihood(temp_chain[1:2,])
      }
      
    } else if (proposal[1,1] >= 0.0 & proposal[2,1] < 0.0){
      proposal[2,1] = temp_chain[2,1]
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[1:2,])
      if (log(runif(1)) < probab){
        temp_chain[1:2,2] = proposal[1:2,1]
        temp_chain[3,2] = bg_likelihood(proposal)
      } else{
        temp_chain[1:2,2] = temp_chain[1:2,1]
        temp_chain[3,2] = bg_likelihood(temp_chain[1:2,])
      }
      
    } else{
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[1:2,])
      if (log(runif(1)) < probab){
        temp_chain[1:2,2] = proposal[1:2,1]
        temp_chain[3,2] = bg_likelihood(proposal)
      } else{
        temp_chain[1:2,2] = temp_chain[1:2,1]
        temp_chain[3,2] = bg_likelihood(temp_chain[1:2,])
      }
    }
    
    # Save the iteration if it is divisible by the divisor without residuals
    if (i%%divisor == 0) {
      chain[,(i/divisor)] = temp_chain[,2]
    }    
    
    # Re-set temporary chain for next iteration
    temp_chain[,1] = temp_chain[,2]
    
    # Print every nth iteration, to know how far along run is
    if (i%%(iterations/20) == 0 | i == 1) {
      print(c(i, "server 6b"))
    }
  }
  return(chain)
}

# Where to start the chain for beta, gamma, and I
startvalue <- array(dim = c(2))
startvalue[1] <- beta # beta guess
startvalue[2] <- gamma  # gamma guess

# Run the MCMC
set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)
rownames(chain) <- c("beta", "gamma", "loglik")

########################################################################################################################

########################################################################################################################

#################
## Saving data ##
#################

beta_gamma_loglik <- t(chain[,]) # transpose beta, gamma, and log likelihood data so that rows become columns etc.
colnames(beta_gamma_loglik) <- c("beta", "gamma", "loglik") # name the columns

# setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")

# Beta, gamma, and log likelihood data
write.csv(data.frame(beta_gamma_loglik), file = "det_mcmc_pop50_R6_g0.15_2_loglik.csv", row.names = FALSE)
