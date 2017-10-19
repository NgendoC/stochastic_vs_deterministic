## Making a stochastic SIS model from scratch
## 11/10/17

library("adaptivetau") # package for stochastic processes

## Set up transitions in and out of each compartment
SI_transitions <- list(
  c(S = -1, I = +1), # infection
  c(I = -1, S = +1) # recovery
)

SI_rateFunc <- function(x, parameters, t) {
  
  # define model parameters in term of the natural parameters
  beta <- parameters["R0"]/parameters["D_inf"] 
  gamma <- 1/parameters["D_inf"]
  
  # create temporary variables for states to simplify the writing of the rates below
  S <- x["S"] 
  I <- x["I"] 
  N <- S + I
  
  # return rates
  return(c(
    beta * S * I / N, # infection 
    gamma * I # recovery
  )) 
}

#############################
### THINGS YOU CAN CHANGE ###
#############################

## Initial values in each compartment
init.values = c(
  S = 10e3-1,
  I = 1
)

## R0 = basic reproduction number, D_inf = duration of infection
parameters <- c(
  R0 = 2, 
  D_inf = 5
)

## Timeframe
tf = 100

par(mfrow=c(2,3))

run_stoch <- ssa.adaptivetau(init.values = init.values, transitions = SI_transitions, 
                       rateFunc = SI_rateFunc, params = parameters, tf = tf)

run_stoch <- data.frame(run)

plot(x = run_stoch$time, y = run_stoch$I, type = "line", col = "red", ylim = c(0,10e3),
  xlab = "Time", ylab = "Number susceptible or infected", main = "Stochastic SIS Model")
  par(new=T)
  plot(x = run_stoch$time, y = run_stoch$S, type = "line", ylim = c(0,10e3), ylab = "", xlab = "") # add susceptible line

## Add legend
legend(80, 1000, c("Susceptible", "Infected"), pch = 1, col = c("black", "red"), bty = "n")