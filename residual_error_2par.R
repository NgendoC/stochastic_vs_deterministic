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
  1, # 5e-3, # beta
  1# 8e-2 # gamma
)

# Optimisation function
sse_fit <- optim(param, sse)
sse_point_data = array((c(sse_fit$par[1], sse_fit$par[2], sse_fit$value)), dim = (c(1,length(param)+1)))
colnames(sse_point_data) <- c("beta", "gamma", "RE")
rownames(sse_point_data) <- c("point")

############################
## Bootstrap optimisation ##
############################

iterations = 20 # how many runs for bootstraps

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

###########
## Plots ##
###########

# Histogram
par(mfrow = c(1,2))

hist(sse_bootstrap_data[,1],nclass=30, main="Beta", xlab="Beta value")
abline(v = sse_fit$par[1], col = "red")

hist(sse_bootstrap_data[,2],nclass=30, main="Gamma", xlab="Gamma value")
abline(v = sse_fit$par[2], col = "red")


# Lines
run_det <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = sse_fit$par))

par(mfrow = c(1,1))
plot(run_stoch$R, ylim = c(0, N), type = "l", col = "orange", xlab = "Timestep", ylab = "Number of individuals")
lines(run_det$I, type = "l", col = "red", xlab = " ", ylab = " ")
lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")
lines(run_det$R, type = "l", col = "black", xlab = "", ylab = "")
legend(100, 0.5*N, c("Deterministic recovered", "True recovered", "Deterministic infected", "True infected"), pch = 1, col = c("black", "orange", "red", "grey"), bty = "n")

# Beta vs. Gamma
par(mfrow = c(1,1))
plot(x = sse_bootstrap_data[,2], y = sse_bootstrap_data[,1], xlab = "Gamma", ylab = "Beta", pch = 20, cex = 0.8)

sse_data <- rbind(sse_point_data, sse_bootstrap_data)
########################################################################################################################

########################################################################################################################

setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")

# Beta, gamma, and residual error data
write.csv(data.frame(sse_data), file = "re_betagamma_test.csv", row.names = FALSE) # point estimate, bootstrap and residual error data

check_data <- read.csv("re_betagamma_test.csv")