## Making a deterministic SIS model from scratch
## 11/10/17

library("deSolve") #package for solving differential equations

## Make an SIS function
si <- function(time, state, parameters) {
  
  # define model parameters in term of the natural parameters
  beta <- parameters["R0"]/parameters["D_inf"] 
  gamma <- 1/parameters["D_inf"]
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I + gamma * I
    dI <-  beta * S * I - gamma * I
    
    return(list(c(dS, dI)))
  })
}

#############################
### THINGS YOU CAN CHANGE ###
#############################

# Proportion in each compartment at the start
init  <- c(
  S = 1-1e-6, 
  I = 1e-6
) # N = 1

## R0 = basic reproduction number, D_inf = duration of infection
parameters <- c(
  R0 = 10, 
  D_inf = 2
)

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
legend(80, 1.0, c("Susceptible", "Infected"), pch = 1, col = c("black","red"), bty = "n")