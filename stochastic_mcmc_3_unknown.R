# MCMC estimating stochastic processes
# 3 unknowns: beta, gamma, and the time of infection for each individual
# 17/11/17

##########################
## Input values for SIR ##
##########################

# Time
timestep <- 0.5
end <- 80
times <- seq(0, end, by = timestep)

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
1/gamma

###############
## The model ##
###############

# Array for holding collective disease status information for whole period of time
data <- array(0, dim =c(length(times), length(init.values)+3))
data[,1] <- times # make first column the timesteps to make plotting easier later on

set.seed(7) # Good ones: 7, 14, 22, 30

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
plot(x = run_stoch$time, y = run_stoch$I, type = "l", col = "red", ylim = c(0,N),
     xlab = "Time", ylab = "Number susceptible/infected/recovered", main = "Stochastic SIR Model")
par(new=T)
plot(x = run_stoch$time, y = run_stoch$S, type = "l", ylim = c(0,N), ylab = "", xlab = "") # add susceptible line
par(new=T)
plot(x = run_stoch$time, y = run_stoch$R, type = "l", col = "orange", ylim = c(0,N), ylab = "", xlab = "") # recovered

# Add legend
legend(60, 0.8*N, c("Susceptible", "Infected", "Recovered"), pch = 1, col = c("black", "red", "orange"), bty = "n")

#################
## Guess new_I ##
#################

inf_period <- ceiling(1/gamma) # mean infectious period calculated from gamma
inf_timestep <- inf_period/timestep # translates infectious days into no. of timesteps

times = times + inf_timestep
zero_array <- array(0, dim = c(inf_timestep, ncol(run_stoch)))
colnames(zero_array) <- c("time","S", "I", "R", "new_I", "new_R")
run_stoch <- rbind(zero_array, run_stoch)
run_stoch$time <- seq(0, (end+inf_period), by = timestep)

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

plot(run_stoch$time, run_stoch$R, ylim = c(0,N), type = "l", col = "orange", xlab = "time (days)", ylab = "Number infectious/recovered")
par(new=T)
plot(run_stoch$time, run_stoch$guess_I, ylim = c(0,N), type = "l", col = "red", xlab = " ", ylab = " ")
par(new=T)
plot(x = run_stoch$time, y = run_stoch$I, type = "l", col = "black", ylim = c(0,N), xlab = " ", ylab = " ")
legend(60, 0.8*N, c("Recovered", "Guessed infected", "True infected"), pch = 1, col = c("orange", "red", "black"), bty = "n")

##################################
## Likelihood, prior, posterior ##
##################################

# Likelihood distribution
likelihood <- function(param){
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  I = as.numeric(param[,2]) # Infected
  new_I = as.numeric(param[,3]) # Newly infected

  total = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:nrow(run_stoch)){
    S = (N - I[i] - run_stoch$R[i]) # Susceptibles
    betalikelihood = dbinom(new_I[i+1], S, (1-(exp(-beta*I[i]*timestep))), log = T)
    gammalikelihood = dbinom(run_stoch$new_R[i+1], I[i], (1-(exp(-gamma*timestep))), log = T)
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
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  
  betaprior = dunif(beta, min = 0, max = 100, log = T)
  gammaprior = dunif(gamma, min = 0, max = 100, log = T)
  return(betaprior + gammaprior)
}

# Posterior distribution
posterior <- function(param){
  return (likelihood(param) + prior(param))
}

# Proposal function for infectious
inf_proposalfunction <- function(param){
  changed_I <- sample(nrow(run_stoch), 1)
  inf_list <- c(-1, 1) # used for choosing -1 or +1 randomly
  
  inf <- sample(c(-1, 1), 1) # will the change at that timepoint be + or - 1 I
  
  neighbour <- if (changed_I == 1){
    changed_I + 1
    } else if (changed_I == nrow(run_stoch)){
      changed_I - 1
      } else{
        changed_I + sample(inf_list, 1) # will choose which neighbouring timepoint is also affected
        }

  param[changed_I[1],2] = param[changed_I[1], 2] + inf
  param[changed_I[1],3] = param[changed_I[1], 3] + inf
  param[neighbour,2] = param[neighbour, 2] - inf
  param[neighbour,3] = param[neighbour, 3] - inf
  return(param)
}

# Proposal function for beta and gamma
proposalfunction <- function(param){
  param[1,1] = rnorm(1, mean = param[1,1], sd = 0.005 ) # beta proposal
  param[2,1] = rnorm(1, mean = param[2,1], sd = 0.1) # gamma proposal
  return(param)
  }

##########
## MCMC ##
##########

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(nrow(startvalue), iterations+1, ncol(startvalue))) # Array for storing chain data

  chain[1,1,1] = startvalue[1,1] # beta
  chain[2,1,1] = startvalue[2,1] # gamma
  chain[,1,2] = startvalue[,2] # all infectious
  chain[,1,3] = startvalue[,3] # newly infectious
  
  for (i in 1:iterations){
    
    # Infectious
    inf_proposal = inf_proposalfunction(chain[,i,])

    # Check susceptibles for negative numbers
    S = array(dim = c(nrow(startvalue)))

    for (time in (1:nrow(startvalue))){
      S[time] = (N - (inf_proposal[time,2] + run_stoch$R[time]))
    }
    
    if(min(inf_proposal[,2]) < 0 | min(inf_proposal[,3]) < 0 | min(S) < 0){
      chain[,i+1,2:3] = chain[,i,2:3]

      } else{
        inf_probab = posterior(inf_proposal) - posterior(chain[,i,])
          if (log(runif(1)) < inf_probab){
            chain[,i+1,2:3] = inf_proposal[,2:3]
            } else{
              chain[,i+1,2:3] = chain[,i,2:3]
              }
        }
    
    # Beta and gamma
    proposal = proposalfunction(chain[,i,]) # beta proposal and gamma proposal
    
    if (proposal[1,1] < 0.0 & proposal[2,1] < 0.0 ){
      chain[1:2,i+1,1] = chain[1:2,i,1]
      
      } else if (proposal[1,1] < 0.0 & proposal[2,1] >= 0.0){
        proposal[1,1] = chain[1,i,1]
        probab = posterior(proposal) - posterior(chain[,i,])
        if (log(runif(1)) < probab){
          chain[,i+1,1] = proposal[,1]
          } else{
            chain[,i+1,1] = chain[,i,1]
            }
      
        } else if (proposal[1,1] >= 0.0 & proposal[2,1] < 0.0){
          proposal[2,1] = chain[2,i,1]
          probab = posterior(proposal) - posterior(chain[,i,])
          if (log(runif(1)) < probab){
            chain[,i+1,1] = proposal[,1]
          } else{
            chain[,i+1,1] = chain[,i,1]
            }
      
          } else{
            probab = posterior(proposal) - posterior(chain[,i,])
            if (log(runif(1)) < probab){
              chain[,i+1,1] = proposal[,1]
              } else{
                chain[,i+1,1] = chain[,i,1]
                }
            }
    # Print every nth iteration, to know how far along run is
    if (i%%(iterations/20) == 0) {
      print(i)
    }
    
  }
  return(chain)
}

# Where to start the chain for beta, gamma, and I
startvalue <- array(dim = c(nrow(run_stoch), 3))
startvalue[1,1] <- 0.01 # beta guess
startvalue[2,1] <- 0.01 # gamma guess
startvalue[,2] <- run_stoch$guess_I # I guess
startvalue[,3] <- run_stoch$guess_new_I # new I guess

# Number of runs
iterations = 10000

# Run the MCMC
set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.5*iterations
acceptance <- 1-mean(duplicated(chain[,-(1:burnIn),]))
inf_acceptance <- 1-mean(duplicated(chain[,-(1:burnIn),2]))

################
## MCMC Plots ##
################

par(mfrow = c(2,2))

hist(chain[1,-(1:burnIn),1],nclass=30, main="Posterior of beta")
abline(v = mean(chain[1,-(1:burnIn),1]), col = "red")

hist(chain[2, -(1:burnIn),1],nclass=30, main="Posterior of gamma")
abline(v = mean(chain[2,-(1:burnIn),1]), col = "red")

plot(chain[1, -(1:burnIn),1], type = "l", main = "Chain values of beta")

plot(chain[2, -(1:burnIn),1], type = "l", main = "Chain values of gamma")

# Plot beta vs. gamma
par(mfrow = c(1,1))
library(RColorBrewer)
library(MASS)

# plot(x = chain[2,,1], y = chain[1,,1], xlab = "Gamma", ylab = "Beta", pch = 20, cex = 0.8)

k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
z <- kde2d(chain[2,,1], chain[1,,1], n=50)
filled.contour(z, nlevels=k, col=my.cols, xlab = "Gamma", ylab = "Beta")
 
par(mfrow = c(1,1))
  
plot(chain[,1,2], ylim = c(0, N), type = "l", col = "red", xlab = "Timestep", ylab = "Number of individuals infected")
  lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")  
  lines(chain[,iterations,2], type = "l", lty = 2, col = "black", xlab = " ", ylab = " ")
legend(130, 1.0*N, c("True infected", "Guessed infected", "MCMC"), pch = 1, col = c("grey", "red", "black"), bty = "n")
  

########################################################################################################################

########################################################################################################################