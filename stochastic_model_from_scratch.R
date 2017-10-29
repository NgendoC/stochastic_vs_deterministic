## Individual-based stochastic SIR model from scratch
## 26/10/17

##################
## Input values ##
##################

# Time
timestep <- 1
#times <- seq(0, 100, by = timestep)
times <- list(0:100)

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

# Array for holding people's disease statuses within a timepoint
status <- array(0, dim = c(N,1))

# Array for holding collective disease status information for whole period of time
data <- array(0, dim =c(length(times[[1]]), length(init.values)))

# For loops for calculating each individual's disease status at each timepoint
for (time in times){
  if (time == 0){ # Set up the number of infecteds and susceptibles at time 0
    status[1:init.values["I"]] <- "I"
    status[init.values["I"]+1:N] <- "S"
  } else{
    for (ind in 1:N){

    # If individual's status at ind is R, they stay R
    if (status[ind] == "R"){ 
      status[ind] <- "R" # they stay recovered
    }
      
    # If status at ind is I:
    } else if (status[ind] == "I"){
    
    r_t <- (1/D_inf)*timestep # Recovery rate
    x <- rbinom(1, N, r_t)
      if (x < r_t){
        status[ind] <- "R"
      } else {
        status[ind] <- "I"
      }
    } else { # if status[ind] == "S"
      # Force of infection at time t
      p <- R0 * (timestep/(D_inf*init.values["N"])) # probability of effective contact
      foi_t <- 1 - (1 - p)^data[time-1, 2] # take data on the infectious from prev timestep!
      # risk of individual becoming infected in next time interval
      x <- runif(1, 0.0, 1.0) # choose a random number between 0 and 1
      if (x < foi_t){ # if x < foi_t, individual becomes infected
        status[ind] <- "I"
      } else { # else individual remains susceptible
        status[ind] <- "S"
      }
  }}
  # Sum up  S/I/R from the status array and add them on into the data array
  data[time, 1] <- length(status[status=="S"])
  data[time, 2] <- length(status[status=="I"])
  data[time, 3] <- length(status[status=="R"])
}
