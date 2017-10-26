## Stochastic SIR model from scratch
## 26/10/17

library("deSolve") #package for solving differential equations

## Make an SIR function
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
  S = 1000-1, 
  I = 1,
  R = 0
) # N = 1

# Make a function for drawing new parameters every time
random_R0 <- function(time){ 
  R0 = rnorm(1, 2, 0.5) # R0
  return(as.numeric(R0))
} 

random_dinf <- function(time){
  dinf = rnorm(1, 5, 0.5) # Duration of infectiousness
  return(as.numeric(dinf))
}

## R0 = basic reproduction number, D_inf = duration of infection
parameters <- c(
  R0 = random_R0,
  D_inf = random_dinf
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
        xlab = "Time", ylab = "Number susceptible or infected",  main = "Stochastic SIR Model",
        lwd = 1, lty = 1, bty = "l", col = c("black","red","orange"))

## Add legend
legend(80, 800, c("Susceptible", "Infected", "Recovered"), pch = 1, col = c("black","red","orange"), bty = "n")