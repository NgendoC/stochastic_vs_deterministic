## Making a deterministic SIS model from scratch
## 11/10/17

library("deSolve") #package for solving differential equations

## Make an SIS function
sir <- function(time, state, parameters) {
  
  # define model parameters in term of the natural parameters
  beta <- parameters["R0"]/parameters["D_inf"] 
  gamma <- 1/parameters["D_inf"]
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I 
    dI <-  beta * S * I - gamma * I
    dR <-  gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

#############################
### THINGS YOU CAN CHANGE ###
#############################

# Proportion in each compartment at the start
init  <- c(
  S = 1-1e-4, 
  I = 1e-4,
  R = 0
) # N = 1

## R0 = basic reproduction number, D_inf = duration of infection
parameters <- c(
  R0 = 2, 
  D_inf = 5
)

## Timeframe
times <- seq(0, 100, by = 1)

## Solve using General Solver for Ordinary Differential Equations (ode)
run <- ode(y = init, times = times, func = sir, parms = parameters)
run_det <- as.data.frame(run) # change to data frame
run_det$time <- NULL # Delete time variable

## Plot model
par(mfrow=c(1,1))

matplot(x = times, y = run_det, type = "l",
  xlab = "Time", ylab = "Proportion susceptible or infected",  main = "Deterministic SIR Model",
  lwd = 1, lty = 1, bty = "l", col = c("black","red","orange"))

## Add legend
legend(80, 1.0, c("Susceptible", "Infected"), pch = 1, col = c("black","red","orange"), bty = "n")