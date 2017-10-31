## Individual-based stochastic SIR model from scratch
## 26/10/17

##################
## Input values ##
##################

# Time
timestep <- 1
times <- seq(0, 100, by = timestep)
#times <- list(0:100)

# Initial population: N-1 susceptible, 1 infectious, 0 recovered
init.values = c(
  S = 10e3-1,
  I = 1,
  R = 0
)
N = sum(init.values)

# R0 & duration of infectiousness
R0 <- 5
D_inf <- 2

# Array for holding collective disease status information for whole period of time
data <- array(0, dim =c(length(times), length(init.values)+1))
data[,1] <- times # make first column the timesteps so you know what's happening when making plots

# Calculating probabilities that do not change with time
p <- R0 * (timestep/(D_inf*N)) # probability of effective contact
r_t <- (1/D_inf)*timestep # Recovery rate

# For loops for calculating each individual's disease status at each timepoint
for (time in times){
  if (time == 0){ # Set up the number of S/I/R at time 0
    data[1,2] <- init.values["S"] # number of susceptibles at time 0
    data[1,3] <- init.values["I"] # number of infecteds at time 0
    data[1,4] <- init.values["R"] # number of recovereds at time 0
    
  } else{
    foi_t <- 1 - (1 - p)^data[time, 3] # take data on the number infectious from previous timestep!
    inf <- rbinom(1, size = data[time,2], prob = foi_t) # number who become infected in this timestep
    rec <- rbinom(1, size = data[time,3], prob = r_t)# number who become recovered in this timestep
    
    data[time+1,2] <- data[time, 2] - inf # number of susceptibles at other times
    
    data[time+1,3] <- data[time, 3]  + inf - rec # number of infecteds at other times
    
    data[time+1,4] <- data[time,4] + rec # number of recovereds at other times
  }
}

run_stoch <- data.frame(data) # make array into a dataframe

plot(x = run_stoch$X1, y = run_stoch$X3, type = "line", col = "red", ylim = c(0,10e3),
     xlab = "Time", ylab = "Number susceptible/infected/recovered", main = "Stochastic SIR Model")
par(new=T)
plot(x = run_stoch$X1, y = run_stoch$X2, type = "line", ylim = c(0,10e3), ylab = "", xlab = "") # add susceptible line
par(new=T)
plot(x = run_stoch$X1, y = run_stoch$X4, type = "line", col = "orange", ylim = c(0,10e3), ylab = "", xlab = "") # recovered

## Add legend
legend(80, 8000, c("Susceptible", "Infected", "Recovered"), pch = 1, col = c("black", "red", "orange"), bty = "n")

######################################################################################################################