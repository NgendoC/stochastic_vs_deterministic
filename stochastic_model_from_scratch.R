## Population-based stochastic SIR model from scratch
## 26/10/17

##################
## Input values ##
##################

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

# R0 & duration of infectiousness
R0 <- 5
D_inf <- 2

# Invent some epidemic data to fit my model to
true_infectious <- (dnorm(times, 20, 4))*N
plot(true_infectious)

#####################
## Values for MCMC ##
#####################

true_mean <- 5
true_var <- 20
true_sd <- sqrt(true_var)
n_obs <- 31

# create values
D <- rnorm(n=n_obs,mean=true_mean,sd=true_sd)

###############
## The model ##
###############

# Array for holding collective disease status information for whole period of time
data <- array(0, dim =c(length(times), length(init.values)+1))
data[,1] <- times # make first column the timesteps to make plotting easier later on

# Calculating probabilities that do not change with time
p <- R0 * (timestep/(D_inf*N)) # probability of effective contact
r_t <- (1/D_inf)*timestep # Recovery rate

# For loops for calculating the numbers susceptible, infected, and recovered at each timepoint
for (time in times){
  if (time == 0){ # Set up the number of S/I/R at time 0
    data[1,2] <- init.values["S"] # number of susceptibles at time 0
    data[1,3] <- init.values["I"] # number of infecteds at time 0
    data[1,4] <- init.values["R"] # number of recovereds at time 0
    
  } else{
    whole_time <- 1/timestep * time # makes time into the whole number that it corresponds to in the array
    foi_t <- 1 - (1 - p)^data[whole_time, 3] # Force of infection at time t, affected by the number infectious at prev t
    inf <- rbinom(1, size = data[whole_time,2], prob = foi_t) # number who become infected in this timestep
    rec <- rbinom(1, size = data[whole_time,3], prob = r_t)# number who become recovered in this timestep
    
    data[whole_time+1,2] <- data[whole_time, 2] - inf # number of susceptibles at other times
    
    data[whole_time+1,3] <- data[whole_time, 3]  + inf - rec # number of infecteds at other times
    
    data[whole_time+1,4] <- data[whole_time,4] + rec # number of recovereds at other times
  }
}

##################################
## Likelihood, prior, posterior ##
##################################

likelihood <- function(param){
  mean = param[1]
  var = param[2]
  
  singlelikelihoods = dnorm(D, mean = mean, sd = sqrt(var), log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

# Prior distribution
prior <- function(param){
  mean = param[1]
  var = param[2]
  meanprior = dnorm(mean, mean = 0, sd = 1, log = T)
  varprior = dnorm(var, mean = 0, sd = 1, log = T)
  return(meanprior+varprior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

##########
## MCMC ##
##########

proposalfunction <- function(param){
  return(rnorm(2,mean = param, sd= c(0.1,0.3)))
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
# Takes a random number from the (Normal) prior distribution
startvalue <- c(1,1) #rnorm(2, prior_mean, prior_sd)

# Number of runs
iterations = 10000
set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.1*iterations
acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

###############
## SIR plots ##
###############

run_stoch <- data.frame(data) # make array into a dataframe
colnames(run_stoch) <- c("time","S", "I", "R")

par(mfrow = c(1,2))

# Plot for SIR model
plot(x = run_stoch$time, y = run_stoch$I, type = "line", col = "red", ylim = c(0,N),
     xlab = "Time", ylab = "Number susceptible/infected/recovered", main = "Stochastic SIR Model")
par(new=T)
plot(x = run_stoch$time, y = run_stoch$S, type = "line", ylim = c(0,N), ylab = "", xlab = "") # add susceptible line
par(new=T)
plot(x = run_stoch$time, y = run_stoch$R, type = "line", col = "orange", ylim = c(0,N), ylab = "", xlab = "") # recovered

# Plot of model fit to data
plot(x = run_stoch$time, y = run_stoch$I, type = "line", col = "red", ylim = c(0,N),
     xlab = "Time", ylab = "Number susceptible/infected/recovered", main = "Model fit to data")
par(new=T)
plot(x = run_stoch$time, y = true_infectious, ylim = c(0,N), ylab = "", xlab = "", type = "line", col="grey") # add susceptible line

## Add legend
legend(60, 0.8*N, c("Susceptible", "Infected", "Recovered"), pch = 1, col = c("black", "red", "orange"), bty = "n")

################
## MCMC Plots ##
################

par(mfrow = c(2,2))

hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of FOI")
abline(v = mean(chain[-(1:burnIn),1]), col = "red")

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of recovery")
abline(v = mean(chain[-(1:burnIn),2]), col = "red")

plot(chain[-(1:burnIn),1], type = "l", main = "Chain values of FOI")

plot(chain[-(1:burnIn),2], type = "l", main = "Chain values of recovery")
