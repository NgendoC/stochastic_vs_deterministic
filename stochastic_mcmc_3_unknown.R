# MCMC estimating stochastic processes
# 3 unknonws: beta, gamma, and the time of infection for each individual
# 17/11/17

##########################
## Input values for SIR ##
##########################

# Time
timestep <- 1
times <- seq(0, 100, by = timestep)

# Initial population: N-1 susceptible, 1 infectious, 0 recovered
init.values = c(
  S = 100-1,
  I = 1,
  R = 0
)
N = sum(init.values)

# Beta & gammma
beta <- 5e-3
gamma <- 8e-2

###############
## The model ##
###############

# Array for holding collective disease status information for whole period of time
data <- array(0, dim =c(length(times), length(init.values)+3))
data[,1] <- times # make first column the timesteps to make plotting easier later on

#set.seed(22) # other good ones: 7, 14, 22, 30

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
    #inf <- rbinom(1, size = data[whole_time,2], prob = foi_t) # number who become infected in this timestep
    #rec <- rbinom(1, size = data[whole_time,3], prob = r_t)# number who become recovered in this timestep
    
    inf <- rbinom(1, size = data[whole_time,2], (1-(exp(-beta*data[whole_time,3]*timestep)))) # number who become infected in this timestep
    rec <- rbinom(1, size = data[whole_time,3], (1-(exp(-gamma*timestep)))) # number who become recovered in this timestep
    
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

par(mfrow = c(1,1))

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
# plot(x = run_stoch$time, y = run_stoch$new_I, type = "line", col = "red", ylim = c(0,max(run_stoch$new_I)),
#      xlab = "Time", ylab = "Number newly infected/recovered", main = "New infections/recoveries")
# par(new=T)
# plot(x = run_stoch$time, y = run_stoch$new_R, type = "line", col = "orange", ylim = c(0,max(run_stoch$new_I)), ylab = "", xlab = "") # add susceptible line

#################
## Guess new_I ##
#################

inf_period <- 2 # infectious period is 2 days

for (i in 1:nrow(run_stoch)){
run_stoch$guess_I[i] <- run_stoch$new_R[i+(inf_period/timestep)] # translate days into timesteps
}

##################################
## Likelihood, prior, posterior ##
##################################

# Likelihood distribution
likelihood <- function(param){
  beta = as.numeric(param[1])
  gamma = as.numeric(param[2])
  # inf = as.numeric(param[3])
  total = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:nrow(run_stoch)){
    betalikelihood = dbinom(run_stoch$new_I[i+1], run_stoch$S[i], (1-(exp(-beta*run_stoch$I[i]*timestep))), log = T)
    gammalikelihood = dbinom(run_stoch$new_R[i+1], run_stoch$I[i], (1-(exp(-gamma*timestep))), log = T)
    
    # betalikelihood = dbinom(run_stoch$new_I[i+1], run_stoch$S[i], (1-(exp(-beta*run_stoch$I[i]*timestep))), log = T)
    # gammalikelihood = dbinom(run_stoch$new_R[i+1], run_stoch$I[i], (1-(exp(-gamma*timestep))), log = T)
    # inflikelihood = dgamma
    total[i] = betalikelihood + gammalikelihood # + inflikelihood
  }
  return(sum(total, na.rm = T))
}

# Prior distribution
prior <- function(param){
  beta = as.numeric(param[1])
  gamma = as.numeric(param[2])
  betaprior = dunif(beta, min = 0, max = 100, log = T)
  gammaprior = dunif(gamma, min = 0, max = 100, log = T)
  return(betaprior + gammaprior)
}

# Posterior distribution
posterior <- function(param){
  return (likelihood(param) + prior(param))
}

proposalfunction <- function(param){ 
  
  beta_prop = rnorm(1, mean = param[1], sd = 0.01)
  gamma_prop = rnorm(1, mean = param[2], sd = 0.01)
  # inf_list <- c(-1, 1)
  # inf <- sample(inf_list, 1)
  
  return(c(beta_prop, gamma_prop)) # , inf))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    if (proposal[1] < 0.0 & proposal[2] < 0.0 ){
      chain[i+1,] = chain[i,]
      
    } else if (proposal[1] < 0.0 & proposal[2] >= 0.0){
      proposal[1] = chain[i,1]
      probab = posterior(proposal) - posterior(chain[i,])
      if (log(runif(1)) < probab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
      }
      
    } else if (proposal[1] >= 0.0 & proposal[2] < 0.0){
      proposal[2] = chain[i,2]
      probab = posterior(proposal) - posterior(chain[i,])
      if (log(runif(1)) < probab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
      }
      
    } else{
      probab = posterior(proposal) - posterior(chain[i,])
      if (log(runif(1)) < probab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
      }
    }
  }
  return(chain)
}

# Where to start the chain
startvalue <- c(0.01,0.01)

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

# Plot beta vs. gamma
par(mfrow = c(1,1))
library(RColorBrewer)
library(MASS)

plot(x = chain[,2], y = chain[,1], xlab = "Gamma", ylab = "Beta", pch = 20, cex = 0.8)
# abline(lm(chain[,1]~chain[,2]), col="red") # regression line

k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
z <- kde2d(chain[,2], chain[,1], n=50)
filled.contour(z, nlevels=k, col=my.cols, xlab = "Gamma", ylab = "Beta")

########################################################################################################################

########################################################################################################################

# A method of randomly choosing -1 or +1 with equal probability
# inf_list <- c(-1, 1)
# sample(inf_list, 1)




