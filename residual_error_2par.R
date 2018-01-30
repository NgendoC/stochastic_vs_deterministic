# Residual error method for comparing a deterministic model with stochastic data
# 2 unknowns: beta, gamma
# 22/01/18

#######################
## Required packages ##
#######################
#install.packages("deSolve")
if (!require("deSolve")) install.packages("deSolve")
library("deSolve") #package for solving differential equations

###############
## Read data ##
###############
setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")
run_stoch <- read.csv("run_stoch.csv")

###########
## Input ##
###########

# Starting point for parameters
beta = 0.005
gamma = 0.08

iterations = 10 # how many runs for bootstraps

#########################
## Deterministic model ##
#########################

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
  
  det_sir <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = param))
  model = det_sir[,4] # Recovered curve for deterministic model
  data = run_stoch[,4] # Recovered curve for stochastic model

  diff_sq = array(0, dim = (c(length(nrow(run_stoch)))))
  
  for (i in 1:nrow(run_stoch)){
    diff_sq[i] = (model[i] - data[i])^2
  }
  
  return(sum(diff_sq))
} 

# Bootstrap function
sse_bootstrap <- function(param){
 
  det_sir <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = param))
  model = det_sir[,4] # Recovered curve for deterministic model
  data = run_stoch[,4] # Recovered curve for stochastic model
  
  data_sample <- sample(1:nrow(run_stoch), nrow(run_stoch), replace=T) # resampling with replacement
  diff_sq = array(0, dim = (c(length(data_sample))))
  
  for (i in data_sample){
    diff_sq[i] = (model[i] - data[i])^2
  }

  return(sum(diff_sq))
} 

##################
## Optimisation ##
##################

# Starting point for parameters
param <- c(
  beta, # beta
  gamma # gamma
)

# Optimisation function
sse_fit <- optim(param, sse)
sse_point_data = array((c(sse_fit$par[1], sse_fit$par[2], sse_fit$value)), dim = (c(1,length(param)+1)))
colnames(sse_point_data) <- c("beta", "gamma", "RE")
rownames(sse_point_data) <- c("point")

############################
## Bootstrap optimisation ##
############################

# Dataset for storing iteration data
sse_bootstrap_data = array(dim = (c(iterations, length(param)+1)))
colnames(sse_bootstrap_data) <- c("beta", "gamma", "RE")

# Bootstrap function for optimisation
for (i in 1:iterations){
  fit <- optim(sse_fit$par, sse_bootstrap) # optimisation function
  sse_bootstrap_data[i,1:2] = fit$par # save beta and gamma values in array
  sse_bootstrap_data[i,3] = fit$value # save residual error
  
  if (i%%(iterations/10) == 0) {
    print(i)
  }
}

########################################################################################################################

########################################################################################################################

#################
## Saving data ##
#################

setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")

# Beta, gamma, and residual error data
sse_data <- rbind(sse_point_data, sse_bootstrap_data)
write.csv(data.frame(sse_data), file = "re_betagamma_test.csv", row.names = FALSE) # point estimate, bootstrap and residual error data
