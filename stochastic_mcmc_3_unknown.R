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
1/gamma

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

#################
## Guess new_I ##
#################

inf_period <- ceiling(1/gamma) # mean infectious period calculated from gamma

for (i in 1:nrow(run_stoch)){
run_stoch$guess_new_I[i] <- run_stoch$new_R[i+(inf_period/timestep)] # guess newly infected, translates days into no. of timesteps
run_stoch$guess_I[i] <- run_stoch$R[i+(inf_period/timestep)] - run_stoch$R[i]
run_stoch$guess_S[i] <- N - (run_stoch$R[i] + run_stoch$guess_I[i])
}
run_stoch[is.na(run_stoch)] <- 0

plot(run_stoch$time, run_stoch$R, ylim = c(0,N), type = "line", col = "orange", xlab = "time (days)", ylab = "Number infectious/recovered")
par(new=T)
plot(run_stoch$time, run_stoch$guess_I, ylim = c(0,N), type = "line", col = "red", xlab = " ", ylab = " ")
par(new=T)
plot(x = run_stoch$time, y = run_stoch$I, type = "line", col = "black", ylim = c(0,N), xlab = " ", ylab = " ")
legend(60, 0.8*N, c("Recovered", "Guessed infected", "True infected"), pch = 1, col = c("orange", "red", "black"), bty = "n")

##################################
## Likelihood, prior, posterior ##
##################################

# Likelihood distribution for beta and gamma
likelihood <- function(param){
  dim(param)
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  I = as.numeric(param[,2])
  new_I = as.numeric(param[,3])
  
  total = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:nrow(run_stoch)){
    betalikelihood = dbinom(new_I[i+1], (N - I[i] - run_stoch$R[i]), (1-(exp(-beta*I[i]*timestep))), log = T)
    gammalikelihood = dbinom(run_stoch$new_R[i+1], I[i], (1-(exp(-gamma*timestep))), log = T)
    total[i] = betalikelihood + gammalikelihood
    # print(betalikelihood)
    # print(gammalikelihood)
  }
  #print(sum(total, na.rm = T))
  return(sum(total, na.rm = T))
}

# Prior distribution
prior <- function(param){
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  
  betaprior = dunif(beta, min = 0, max = 100, log = T)
  gammaprior = dunif(gamma, min = 0, max = 100, log = T)
  return(betaprior + gammaprior)
}

# Posterior distribution for beta and gamma
posterior <- function(param){
  return (likelihood(param) + prior(param))
}

# Proposal function for infectious
inf_proposalfunction <- function(param){
  changed_I <- sample(length(times), 2)
  inf_list <- c(-1, 1) # used for choosing -1 or +1 randomly
  inf <- sample(inf_list, 1) # will the change at that timepoint be + or - 1 I
  
  #param[,] = chain[,i,]
  param[changed_I[1],2] = param[changed_I[1], 2] + inf
  param[changed_I[1],3] = param[changed_I[1], 3] + inf
  param[changed_I[2],2] = param[changed_I[2], 2] - inf # take away infected at previous timestep so altogether infected stay the same
  param[changed_I[2],3] = param[changed_I[2], 3] - inf
  
  return(param)
}

# Proposal function for beta and gamma
proposalfunction <- function(param){
  param[1,1] = rnorm(1, mean = param[1,1], sd = 0.001 ) # beta proposal
  param[2,1] = rnorm(1, mean = param[2,1], sd = 0.01) # gamma proposal
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
    inter_chain = chain[,i,] # intermediate chain for this iteration
    
    inf_proposal = inf_proposalfunction(inter_chain)
    
    if(min(inf_proposal[,2]) < 0 | min(inf_proposal[,3]) < 0){
      chain[,i+1,] = chain[,i,]
        
    } else{ 
      inf_probab = posterior(inf_proposal) - posterior(chain[,i,])
      if (log(runif(1)) < inf_probab){
      chain[,i+1,] = inf_proposal
      }else{
      chain[,i+1,] = chain[,i,]
      }
    }
    
    # Beta and gamma
    proposal = proposalfunction(chain[,i,]) # beta proposal and gamma proposal
    
    if (proposal[1,1] < 0.0 & proposal[2,1] < 0.0 ){
      chain[1:2,i+1,1] = chain[1:2,i,1]
      
    } else if (proposal[1,1] < 0.0 & proposal[2,1] >= 0.0){
      proposal[1,1] = chain[1,i,1]
      probab = posterior(proposal) - posterior(chain[,i,])
      #print(probab)
      if (log(runif(1)) < probab){
        chain[,i+1,] = proposal
      }else{
        chain[,i+1,] = chain[,i,]
      }
      
    } else if (proposal[1,1] >= 0.0 & proposal[2,1] < 0.0){
      proposal[2,1] = chain[2,i,1]
      probab = posterior(proposal) - posterior(chain[,i,])
      #print(probab)
      if (log(runif(1)) < probab){
        chain[,i+1,] = proposal
      }else{
        chain[,i+1,] = chain[,i,]
      }
      
    } else{
      probab = posterior(proposal) - posterior(chain[,i,])
      #print(probab)
      if (log(runif(1)) < probab){
        chain[,i+1,] = proposal
      }else{
        chain[,i+1,] = chain[,i,]
      }
    }
  }
  return(chain)
}

# Where to start the chain for beta, gamma, and I
startvalue <- array(dim = c(length(times), 3))
startvalue[1,1] <- 0.01 # beta guess
startvalue[2,1] <- 0.01 # gamma guess
startvalue[,2] <- run_stoch$guess_I # I guess
startvalue[,3] <- run_stoch$guess_new_I # new I guess

# Number of runs
iterations = 5000

# Run the MCMC
#set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.1*iterations
acceptance <- 1-mean(duplicated(chain[,-(1:burnIn),]))

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

plot(x = chain[2,,1], y = chain[1,,1], xlab = "Gamma", ylab = "Beta", pch = 20, cex = 0.8)
# abline(lm(chain[,1]~chain[,2]), col="red") # regression line

k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
z <- kde2d(chain[2,,1], chain[1,,1], n=50)
filled.contour(z, nlevels=k, col=my.cols, xlab = "Gamma", ylab = "Beta")

# plot(run_stoch$time, run_stoch$R, ylim = c(0, N), type = "line", col = "orange", xlab = "time (days)", ylab = "Number infectious/recovered")
# par(new=T)
# plot(chain[,1,2], ylim = c(0, N), type = "line", col = "red", xlab = " ", ylab = " ")
# par(new=T)
# plot(x = run_stoch$time, y = run_stoch$I, ylim = c(0, N), type = "line", col = "black", xlab = " ", ylab = " ")
# legend(60, 0.8*N, c("Recovered", "Guessed infected", "True infected"), pch = 1, col = c("orange", "red", "black"), bty = "n")
 
plot(chain[,1,2], ylim = c(0, N), type = "line", col = "red", xlab = " ", ylab = " ")
par(new=T)
plot(chain[,5,2], ylim = c(0, N), col = "green", lty = "2", xlab = " ", ylab = " ")
par(new=T)
plot(x = run_stoch$time, y = run_stoch$I, ylim = c(0, N), type = "line", col = "black", xlab = " ", ylab = " ")
legend(60, 0.8*N, c("Recovered", "Guessed infected", "True infected"), pch = 1, col = c("orange", "red", "black"), bty = "n")
########################################################################################################################

########################################################################################################################

# numbers <- array(dim = c(3,3,3))
# numbers[1,,] <- 1:3
# numbers[2,,] <- -1:1
# numbers[3,,] <- 4:6
# 
# min(numbers[1,,])
