## Making a stochastic SIS model from scratch
## 11/10/17

library("adaptivetau") # package for stochastic processes

## Set up initial values in each compartment
init.values = c(
  S = 999,
  I = 1
)

## Set up transitions in and out of each compartment
SI_transitions <- list(
  c(S = -1, I = +1), # infection
  c(I = -1, S = +1) # recovery
)

SI_rateFunc <- function(x, parameters, t) {
  
  # define model parameters
  beta <- 0.5 
  gamma <- 0.2
  
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

parameters <- c(beta = 0.5, gamma = 0.2)
run <- ssa.adaptivetau(init.values = c(S = 999, I = 1), transitions = SI_transitions, 
                       rateFunc = SI_rateFunc, params = parameters, tf = 100)

run_df <- data.frame(run)

plot(x = run_df$time, y = run_df$I, type = "line", col = "red", ylim = c(0,1000),
  xlab = "Time", ylab = "Number susceptible or infected", main = "Stochastic SIS Model")
  par(new=T)
  plot(x = run_df$time, y = run_df$S, type = "line", ylim = c(0,1000), ylab = "") # add susceptible line

## Add legend
legend(80, 1000, c("Susceptible", "Infected"), pch = 1, col = c("black", "red"), bty = "n")