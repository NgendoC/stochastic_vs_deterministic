# Residual error heatmap
# 01/02/18

#######################
## Required packages ##
#######################
#install.packages("deSolve")
if (!require("deSolve")) install.packages("deSolve")
library("deSolve") #package for solving differential equations

if (!require("plotly")) install.packages("plotly")
library("plotly") #package for solving differential equations

###############
## Read data ##
###############
setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")
run_stoch <- read.csv("run_stoch.csv")

######################################################
## Make all possible combinations of beta and gamma ##
######################################################

# Starting point for parameters
beta <- seq(0.003, 0.008, by = ((0.008-0.003)/50))
gamma <- seq(0.05, 0.20, by = ((0.20-0.05)/50))

bg_combo <- expand.grid(beta, gamma)

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
  
  diff_sq = array(dim = (c(length(nrow(run_stoch)))))
  
  for (i in 1:nrow(run_stoch)){
    diff_sq[i] = (model[i] - data[i])^2
  }
  
  return(sum(diff_sq))
} 

# Dataset for storing data
sse_heatmap_data = array(dim = (c(nrow(bg_combo), 3)))
colnames(sse_heatmap_data) <- c("beta", "gamma", "RE")
sse_heatmap_data <- as.data.frame(sse_heatmap_data)

# Calculating residual errors
for (i in 1:nrow(bg_combo)){
  re <- sse(bg_combo[i, ]) # Gives you the RE for a given combination
  # fit <- optim(param, sse_bootstrap, constant = i) # optimisation function, original guess as starting point
  sse_heatmap_data[i,1:2] = bg_combo[i, ] # save beta and gamma values in array
  sse_heatmap_data[i,3] = re # save residual error
  
  # bootstrap_param[3] = bootstrap_param[3]+1
  if (i%%(nrow(bg_combo)/10) == 0) {
    print(i)
  }
}

#################
## Saving data ##
#################

setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")

# Beta, gamma, and residual error data
write.csv(sse_heatmap_data, file = "re_heatmap_test.csv", row.names = FALSE) # point estimate, bootstrap and residual error data
