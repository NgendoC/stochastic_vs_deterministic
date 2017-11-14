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
  #sum_betagamma = array(dim = c(nrow(run_stoch), 1))
  total = 0
  
  for (i in 1:nrow(run_stoch)){
  betalikelihood = dbinom(run_stoch$new_I, run_stoch$S, (1-(exp(-beta*run_stoch$I*timestep))), log = T)
  gammalikelihood = dbinom(run_stoch$new_R, run_stoch$I, (1-(exp(-gamma*timestep))), log = T)
  #sum_betagamma[i] = betalikelihood + gammalikelihood
  #print(sum_betagamma[i])
  total = total + (betalikelihood + gammalikelihood)
  }
  #total_sum = sum(sum_betagamma)
  #print(total_sum)
  print(sum(total))
  return(sum(total))
}

# sum_betagamma = array(dim = c(nrow(run_stoch), 1))
# 
# for (i in 1:nrow(run_stoch)){
#   betalikelihood = dbinom(run_stoch$new_I, run_stoch$S, (1-(exp(-0.1*run_stoch$I*timestep))), log = T)
#   gammalikelihood = dbinom(run_stoch$new_R, run_stoch$I, (1-(exp(-0.2*timestep))), log = T)
#   sum_ll = betalikelihood + gammalikelihood
#   sum_betagamma[i] = sum_ll
#   print(sum_betagamma[i])
# }

# Prior distribution
prior <- function(param){
  beta = as.numeric(param[1])
  gamma = as.numeric(param[2])
  
  betaprior = dunif(beta, min = 0, max = 100, log = T)
  gammaprior = dunif(gamma, min = 0, max = 100, log = T)
  print(betaprior + gammaprior)
  return(betaprior + gammaprior)
}

# Posterior distribution
posterior <- function(param){
  return (likelihood(param) + prior(param))
}

proposalfunction <- function(param){ # beta and gamma need to be >0
  
  # beta_prop = rnorm(1, mean = param[1], sd = 1)
  # gamma_prop = rnorm(1, mean = param[2], sd = 1)
  beta_prop = runif(1, min = 0, max = 10)
  gamma_prop = runif(1, min = 0, max = 10)
  
  #print(c(beta_prop, gamma_prop))
  return(c(beta_prop, gamma_prop))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    # if (proposal[1] < 0.0 | proposal[2] < 0.0){
    #   chain[i+1,] = chain[i,]
    # }
    # else{
      probab = exp(posterior(proposal) - posterior(chain[i,]))
      if (runif(1) < probab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
        }
      #}
  }
  return(chain)
}

# Where to start the chain
startvalue <- c(0.1,0.1)

# Number of runs
iterations = 10000
#set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.1*iterations
acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

################
## MCMC Plots ##
################

par(mfrow = c(2,2))

hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of beta")
abline(v = mean(chain[-(1:burnIn),1]), col = "red")

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of gamma")
abline(v = mean(chain[-(1:burnIn),2]), col = "red")

plot(chain[-(1:burnIn),1], type = "l", main = "Chain values of beta")

plot(chain[-(1:burnIn),2], type = "l", main = "Chain values of gamma")

########################################################################################################################

########################################################################################################################