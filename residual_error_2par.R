# Residual error method for comparing a deterministic model with stochastic data
# 2 unknowns: beta, gamma
# 22/01/18

#######################
## Required packages ##
#######################
#install.packages("deSolve")
library("deSolve") #package for solving differential equations

##########################
## Input values for SIR ##
##########################

# Time
timestep <- 0.5
end <- 80
times <- seq(0, end, by = timestep)

# Initial population: N-1 susceptible, 1 infectious, 0 recovered
init.values = c(
  S = (100-1),
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

#########################
## Deterministic model ##
#########################

sir <- function(time, state, param) {
  
  # define model parameters in term of the natural parameters
  beta <- param[1] 
  gamma <- param[2]
  
  with(as.list(c(state, param)), {
    
    dS <- -beta * S * I 
    dI <- (beta * S * I) -(gamma * I)
    dR <-  gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

####################
## Residual Error ##
####################

# Function for calculation sum of squared differences (deterministic model vs. stochastic data)
sse <- function(param){
  
  data = run_stoch[,4]
  det_sir <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = param))
  model = det_sir[,4]
  
  diff_sq = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:nrow(run_stoch)){
    diff_sq[i] = (model[i] - data[i])^2
  }
  
  return(sum(diff_sq))
} 

# Starting point for parameters
param <- c(
  5e-2, # 5e-3, # beta
  8e-1# 8e-2 # gamma
)

# Function for optimisation
fit <- optim(param, sse)

# Plot the results
run_det <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = fit$par))

plot(run_stoch$R, ylim = c(0, N), type = "l", col = "orange", xlab = "Timestep", ylab = "Number of individuals")
lines(run_det$I, type = "l", col = "red", xlab = " ", ylab = " ")
lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")
lines(run_det$R, type = "l", col = "black", xlab = "", ylab = "")
legend(100, 0.5*N, c("Deterministic recovered", "True recovered", "Deterministic infected", "True infected"), pch = 1, col = c("black", "orange", "red", "grey"), bty = "n")


# Bootstrapping for confidence intervals

# https://ms.mcmaster.ca/~bolker/eeid/2011_eco/EEID2011_Fitting.pdf