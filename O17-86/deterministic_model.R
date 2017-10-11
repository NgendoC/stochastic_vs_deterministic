## Making a deterministic SIS model from scratch
## 11/10/17

library("deSolve") #package for solving differential equations

## Make an SIS function
si <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I + gamma * I
    dI <-  beta * S * I - gamma * I
    
    return(list(c(dS, dI)))
  })
}

## Set parameters for SI function

# Proportion in each compartment at the start
init  <- c(S = 1-1e-6, I = 1e-6) # N = 1

## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 0.5, gamma = 0.2)

## Timeframe
times <- seq(0, 100, by = 1)

## Solve using General Solver for Ordinary Differential Equations (ode)
out <- ode(y = init, times = times, func = si, parms = parameters)
out <- as.data.frame(out) # change to data frame
out$time <- NULL # Delete time variable

## Plot model
matplot(x = times, y = out, type = "l",
  xlab = "Time", ylab = "Proportion susceptible or infected", main = "Deterministic SIS Model",
  lwd = 1, lty = 1, bty = "l", col = c("black","red"))

## Add legend
legend(60, 1.0, c("Susceptible", "Infected"), pch = 1, col = c("black","red"), bty = "n")

