## Population-based stochastic SIR model from scratch
## 26/10/17

##########################
## Input values for SIR ##
##########################

# Time
timestep <- 0.5
end <- 80
times <- seq(0, end, by = timestep)

# Initial population: N-1 susceptible, 1 infectious, 0 recovered
N <- 1000
start_I <- 20

init.values = c(
  S = N-start_I,
  I = start_I,
  R = 0
)

# Beta & gammma
R0 <- 1.5*4 # R0 = beta*N/gamma
gamma <- 0.15
beta <- R0*gamma/N
print(beta)

###############
## The model ##
###############

# Array for holding collective disease status information for whole period of time
data <- array(0, dim =c(length(times), length(init.values)+3))
data[,1] <- times # make first column the timesteps to make plotting easier later on

# for (i in 1:10000){
  set.seed(3117)
  # R0 = 1.5 || Beta = 0.15, gamma = 0.1
# Pop = 50: 7, 12, 13, 14, 18
# Pop = 200: 14, 16, 18, 22, 30
# Pop = 1000: 3, 14, 30, 35, 41

# R0 = 7.5 || Beta = 0.75, gamma = 0.1
# Pop = 50: 1, 4, 7, 10, 15
# Pop = 200: 2, 7, 13, 15, 16
# Pop = 1000: 1, 3, 7, 15, 17

# Beta = 0.005, gamma = 0.08
# Pop = 50: 35, 57, 95, 113, 227
# Pop = 200: 2, 7, 8, 14, 20
# Pop = 1000: 6, 8, 20, 27, 28

# Beta = 0.05, gamma = 0.08
# Pop = 50: 7, 10, 15, 18, 25
# Pop = 200: 4, 7, 9, 15, 19
# Pop = 1000: 6, 8, 27, 28, 30

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
    
    inf <- rbinom(1, size = data[whole_time,2], (1-(exp((-beta)*data[whole_time,3]*timestep)))) # number who become infected in this timestep
    rec <- rbinom(1, size = data[whole_time,3], (1-(exp((-gamma)*timestep)))) # number who become recovered in this timestep
    
    data[whole_time+1,2] <- data[whole_time,2] - inf # number of susceptibles at next timestep
    
    data[whole_time+1,3] <- data[whole_time,3]  + inf - rec # number of infecteds at next timestep
    
    data[whole_time+1,4] <- data[whole_time,4] + rec # number of recovereds at next timestep
    
    data[whole_time+1,5] <- data[whole_time+1,3] - data[whole_time,3] + data[whole_time+1,4] - data[whole_time,4] # number of newly infected
    
    data[whole_time+1,6] <- data[whole_time+1,4] - data[whole_time,4] # number of newly recovered
  }
}
  # if (data[nrow(data), 3] == 0 & max(data[, 3]) > 0.1*N){
  #   print(i)
  # }
  # }

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

###############
## Save data ##
###############

setwd("/home/evelina/Development/stochastic_vs_deterministic")

# Save SIR data
write.csv(run_stoch, file = "run_stoch.csv", row.names = FALSE)