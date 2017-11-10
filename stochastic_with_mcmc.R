## Stochastic model with MCMC
## 04/11/17

##########################
## Input values for SIR ##
##########################

# Time
timestep <- 0.1
times <- seq(0, 100, by = timestep)

# Initial population: N-1 susceptible, 1 infectious, 0 recovered
init.values = c(
  S = 100-1,
  I = 1,
  R = 0
)
N = sum(init.values)

# Invent some epidemic data to fit my model to
true_infectious <- (dnorm(times, 20, 4))*N
#plot(true_infectious)

###############
## The model ##
###############

# Array for holding collective disease status information for whole period of time
data <- array(0, dim =c(length(times), length(init.values)+3))
data[,1] <- times # make first column the timesteps to make plotting easier later on
  
# R0 & duration of infectiousness
R0 <- 3
D_inf <- 2
  
# Calculating probabilities that do not change with time
p <- R0 * (timestep/(D_inf*N)) # probability of effective contact
r_t <- (1/D_inf)*timestep # Recovery rate
  
# For loops for calculating the numbers susceptible, infected, and recovered at each timepoint
for (time in times){
  if (time == 0){ # Set up the number of S/I/R at time 0
    data[1,2] <- init.values["S"] # number of susceptibles at time 0
    data[1,3] <- init.values["I"] # number of infecteds at time 0
    data[1,4] <- init.values["R"] # number of recovereds at time 0
    data[1,5] <- init.values["I"] # number newly infected at time 0
    data[1,6] <- init.values["R"] # number newly recovered at time 0
      
  } else{
    whole_time <- 1/timestep * time # makes time into the whole number that it corresponds to in the array
    foi_t <- 1 - (1 - p)^data[whole_time, 3] # Force of infection at time t, affected by the number infectious at prev t
    inf <- rbinom(1, size = data[whole_time,2], prob = foi_t) # number who become infected in this timestep
    rec <- rbinom(1, size = data[whole_time,3], prob = r_t)# number who become recovered in this timestep
      
    data[whole_time+1,2] <- data[whole_time,2] - inf # number of susceptibles at other times
      
    data[whole_time+1,3] <- data[whole_time,3]  + inf - rec # number of infecteds at other times
      
    data[whole_time+1,4] <- data[whole_time,4] + rec # number of recovereds at other times
    
    data[whole_time+1,5] <- data[whole_time+1,3] - data[whole_time,3] + data[whole_time+1,4] - data[whole_time,4] # number of newly infected
    
    data[whole_time+1,6] <- data[whole_time+1,4] - data[whole_time,4] # number of newly recovered
  }
}

###############
## SIR plots ##
###############

run_stoch <- data.frame(data) # make array into a dataframe
colnames(run_stoch) <- c("time","S", "I", "R", "new_I", "new_R")

par(mfrow = c(1,2))

# Plot for SIR model
plot(x = run_stoch$time, y = run_stoch$I, type = "line", col = "red", ylim = c(0,N),
     xlab = "Time", ylab = "Number susceptible/infected/recovered", main = "Stochastic SIR Model")
par(new=T)
plot(x = run_stoch$time, y = run_stoch$S, type = "line", ylim = c(0,N), ylab = "", xlab = "") # add susceptible line
par(new=T)
plot(x = run_stoch$time, y = run_stoch$R, type = "line", col = "orange", ylim = c(0,N), ylab = "", xlab = "") # recovered

# Add legend
legend(60, 0.8*N, c("Susceptible", "Infected", "Recovered"), pch = 1, col = c("black", "red", "orange"), bty = "n")

# Plot for newly infected and newly recovered
plot(x = run_stoch$time, y = run_stoch$new_I, type = "line", col = "red", ylim = c(0,max(run_stoch$new_I)),
     xlab = "Time", ylab = "Number newly infected/recovered", main = "New infections/recoveries")
par(new=T)
plot(x = run_stoch$time, y = run_stoch$new_R, type = "line", col = "orange", ylim = c(0,max(run_stoch$new_I)), ylab = "", xlab = "") # add susceptible line

##################################
## Likelihood, prior, posterior ##
##################################

likelihood <- function(param){
  beta = as.numeric(param[1])
  gamma = as.numeric(param[2])
  sumll = 0 #sum of likelihoods is zero before you start counting
  
  for (time in times){
  betalikelihood = dbinom(run_stoch$new_I, run_stoch$S, (1-(exp(-beta*run_stoch$I*timestep))), log = T)
  gammalikelihood = dbinom(run_stoch$new_R, run_stoch$I, (gamma*timestep), log = T)
  #sumll = sum(betalikelihood + gammalikelihood)
  sumll = sumll + (betalikelihood + gammalikelihood)
  }
  #print(betalikelihood)
  #print(gammalikelihood)
  print(sumll)
  return(sumll)
}

# Prior distribution
prior <- function(param){
  beta = as.numeric(param[1])
  gamma = as.numeric(param[2])
  
  betaprior = dunif(beta, min = 0, max = Inf, log = T)
  gammaprior = dunif(gamma, min = 0, max = Inf, log = T)
  # print(betaprior)
  # print(gammaprior)
  return(betaprior + gammaprior)
}

# Posterior distribution
posterior <- function(param){
  return (likelihood(param) + prior(param))
}

proposalfunction <- function(param){
  return(rnorm(2, mean = param, sd = 1))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])

    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }

  return(chain)
}

# Where to start the chain
startvalue <- c(0,0)

# Number of runs
iterations = 10000
set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.1*iterations
acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

################
## MCMC Plots ##
################

# par(mfrow = c(2,2))
# 
# hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of beta")
# abline(v = mean(chain[-(1:burnIn),1]), col = "red")
# 
# hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of gamma")
# abline(v = mean(chain[-(1:burnIn),2]), col = "red")
# 
# plot(chain[-(1:burnIn),1], type = "l", main = "Chain values of beta")
# 
# plot(chain[-(1:burnIn),2], type = "l", main = "Chain values of gamma")

########################################################################################################################

########################################################################################################################

## The original MCMC

# likelihood <- function(param){
#   mean = param[1]
#   var = param[2]
#   
#   singlelikelihoods = dnorm(D, mean = mean, sd = sqrt(var), log = T)
#   sumll = sum(singlelikelihoods)
#   return(sumll)   
# }
# 
# # Prior distribution
# prior <- function(param){
#   mean = param[1]
#   var = param[2]
#   meanprior = dnorm(mean, mean = 0, sd = 1, log = T)
#   varprior = dnorm(var, mean = 0, sd = 1, log = T)
#   return(meanprior+varprior)
# }
# 
# posterior <- function(param){
#   return (likelihood(param) + prior(param))
# }

# proposalfunction <- function(param){
#   return(rnorm(2,mean = param, sd= c(0.1,0.3)))
# }
# 
# run_metropolis_MCMC <- function(startvalue, iterations){
#   chain = array(dim = c(iterations+1,2))
#   chain[1,] = startvalue
#   for (i in 1:iterations){
#     proposal = proposalfunction(chain[i,])
#     
#     probab = exp(posterior(proposal) - posterior(chain[i,]))
#     if (runif(1) < probab){
#       chain[i+1,] = proposal
#     }else{
#       chain[i+1,] = chain[i,]
#     }
#   }
#   
#   return(chain)
# }


########################
## Rate of infection  ##
########################

# inf_likelihood <- function(inf_param){
#   inf_param = foi_t # force of infection at time t
#   inf_singlell = dbinom(inf_data, size=S(t-1), prob = foi_t, log = T) # inf_param is in the foi_t
#   inf_sumll = sum(singlell)
#   return(inf_sumll)
#}

# inf_prior <- function(inf_param){
#   infprior = dbeta(inf_param, shape1, shape2, ncp = 0, log = T)
#   return(infprior)
#}

# inf_posterior <- function(inf_param){
#   return(inf_likelihood + inf_prior)
#}

####################
## Recovery rate  ##
####################

#  rec_likelihood <- function(rec_param){
#   rec_param = r_t # recovery rate at time t
#   rec_singlell = dbinom(rec_data, size=I(t-1), prob = r_t, log = T) # rec_param is in the r_t
#   rec_sumll = sum(rec_singlell)
#   return(rec_sumll)
#}

# rec_prior <- function(rec_param){
#   recprior = dbeta(rec_param, shape1, shape2, ncp = 0, log = T)
#   return(recprior)
#}

# rec_posterior <- function(rec_param){
#   return(rec_likelihood + rec_prior)
#}

##########
## MCMC ##
##########

# proposalfunction <- function(param){
#   return(rnorm(1,mean = param, sd= 1))
# }

# run_metropolis_MCMC <- function(startvalue, iterations){
#   chain = array(dim = c(iterations+1,2))
#   chain[1,] = startvalue
#   
#   for (i in 1:iterations){
#     # Infection rate
#     inf_proposal = proposalfunction(chain[i,1])
#     inf_probab = exp(inf_posterior(inf_proposal) - inf_posterior(chain[i,1]))
#     if (runif(1) < inf_probab){
#       chain[i+1,1] = inf_proposal
#     }else{
#       chain[i+1,1] = chain[i,1]
#     }
#     # Recovery rate
#     rec_proposal = proposalfunction(chain[i,])
#     rec_probab = exp(rec_posterior(rec_proposal) - rec_posterior(chain[i,2]))
#     if (runif(1) < rec_probab){
#       chain[i+1,2] = rec_proposal
#     }else{
#       chain[i+1,2] = chain[i,2]
#     }
#   }
#   
#   return(chain)
# }