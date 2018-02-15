# Residual error method for comparing a deterministic model with stochastic data
# 2 unknowns: beta, gamma
# 22/01/18

#######################
## Required packages ##
#######################
# install.packages("deSolve")
# if (!require("deSolve")) install.packages("deSolve")
library("deSolve") #package for solving differential equations

###############
## Read data ##
###############
setwd("/home/evelina/Development/stochastic_vs_deterministic")
run_stoch <- read.csv("data_pop1000_b0.75_g0.1_1.csv")

###########
## Input ##
###########

# Starting point for parameters
R0 <-3.0 # R0 = beta*N/gamma
gamma <- 0.1
beta <- R0*gamma/N

iterations = 1000 # how many runs for bootstraps

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
    
    dS <- -(beta * S * I) 
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

# Bootstrap function
sse_bootstrap <- function(param, constant){
  
  det_sir <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = param))
  model = det_sir[,4] # Recovered curve for deterministic model
  data = run_stoch[,4] # Recovered curve for stochastic model
  
  set.seed(constant) # keeps the data_sample number set constant for the iteration
  data_sample <- sample(1:nrow(run_stoch), nrow(run_stoch), replace=T) # resampling with replacement
  diff_sq = array(dim = (c(length(data_sample))))
  
  for (i in 1:length(data_sample)){
    diff_sq[i] = (model[data_sample[i]] - data[data_sample[i]])^2
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
  constant = i
  # fit <- optim(sse_fit$par, sse_bootstrap, constant = i) # optimisation function, point estimate as starting point
  fit <- optim(param, sse_bootstrap, constant = i) # optimisation function, original guess as starting point
  sse_bootstrap_data[i,1:2] = fit$par[1:2] # save beta and gamma values in array
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

setwd("/home/evelina/OneDrive/MRes_BMR/Project_1/Work_folder/Data")

# Beta, gamma, and residual error data
sse_data <- rbind(sse_point_data, sse_bootstrap_data)
# write.csv(data.frame(sse_data), file = "re_pop1000_b0.75_g0.1_test.csv", row.names = FALSE) # point estimate, bootstrap and residual error data

# Plots
par(mfrow = c(1,2))
hist(sse_data[2:nrow(sse_data),1],nclass=30, col="gray", main="RE beta seed 1", xlab="")
abline(v = sse_data[1,1], col = "red")
box()
hist(sse_data[2:nrow(sse_data),2],nclass=30, col="gray", main="RE gamma seed 1", xlab="")
abline(v = sse_data[1,2], col = "red")
box()

par(mfrow = c(1,1))
hist(sse_data[2:nrow(sse_data), 3],nclass=30, col="gray", main="RE seed 1", xlab="Residual Error")
abline(v = sse_data[1,3], col = "red")
box()

par(mfrow= c(1,2))
plot(run_stoch$I, ylim = c(0, N), type = "l", col = "white", xlab = "Timestep", ylab = "Number of individuals infected", main = "Point estimate")
  det_sir <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = sse_data[1,1:2]))
  lines(run_stoch$I, ylim = c(0, 50), type = "l", col = "red", lwd = 3)
  lines(run_stoch$R, ylim = c(0, 50), type = "l", col = "orange", lwd = 3)
  lines(det_sir$I, type = "l", lty = 2, col = "black", xlab = " ", ylab = " ", lwd = 1)
  lines(det_sir$R, type = "l", lty = 2, col = "black", xlab = " ", ylab = " ", lwd = 1)
  
  
  lines(run_stoch$I, ylim = c(0, 50), type = "l", col = "red", lwd = 3)
  legend(100, 0.8*N, c("True infected", "Deterministic infected"), pch = 1, col = c("red", "black"), bty = "n")
  plot(run_stoch$I, ylim = c(0, N), type = "l", col = "white", xlab = "Timestep", ylab = "Number of individuals infected", main = "Bootstrap")
  lines(run_stoch$I, ylim = c(0, 50), type = "l", col = "red", lwd = 3)
  lines(run_stoch$R, ylim = c(0, 50), type = "l", col = "orange", lwd = 3)  
  for (i in 2:nrow(sse_data)){
    det_sir <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = sse_data[i,1:2]))
    lines(det_sir$I, type = "l", lty = 2, col = "black", xlab = " ", ylab = " ", lwd = 0.5)
    lines(det_sir$R, type = "l", lty = 2, col = "black", xlab = " ", ylab = " ", lwd = 0.5)
  }
  legend(100, 0.8*N, c("True infected", "Deterministic infected"), pch = 1, col = c("red", "black"), bty = "n")
  