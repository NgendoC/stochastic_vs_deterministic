## Individual-based stochastic SIR model from scratch
## 26/10/17

# Input values

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

# Time
timestep <- 1
times <- seq(0, 100, by = timestep)

# Array for holding people's disease statuses
status <- array(0, dim = c(N,1))

for time in times{
  for ind in 1:N{
    # Force of infection at time t
    p <- R0 * (timestep/(D_inf*init.values["N"])) # probability of effective contact
    foi_t <- 1 - (1 - p)^dI
     # risk of individual becoming infected in next time interval
    x <- runif(1, 0.0, 1.0) # choose a random number between 0 and 1
    # if x < foi_t, individual becomes infected
    # else individual remains susceptible

    # Recovery 
    r_t <- (1/D_inf)*timestep
    rbinom(1, N, r_t)
  
}}
