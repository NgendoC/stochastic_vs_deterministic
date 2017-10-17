## Deterministic & stochastic SIS model
## 12/10/2017

## Necessary libraries
library("deSolve")       # differential equation solver for deterministic model
library("adaptivetau")   # for stochastic processes

## Initial values in each compartment
stoch_init.values = c(
  S = 10000-1,               # Susceptible
  I = 1                  # Infected
)

det_init.values = c(
  S = 1-10e-5,               # Susceptible
  I = 10e-5                  # Infected
)

######################
## Stochastic model ##
######################

## Set up transitions in and out of each compartment
SIS_transitions <- list(
  c(S = -1, I = +1), # infection
  c(I = -1, S = +1) # recovery
)

stoch_SIS_rateFunc <- function(x, parameters, t) {
  
  # define model parameters
  beta <- 1.0
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

#########################
## Deterministic model ##
#########################

det_SIS_rateFunc <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I + gamma * I
    dI <-  beta * S * I - gamma * I
    
    return(list(c(dS, dI)))
  })
}

parameters <- c(beta = 1.0, gamma = 0.2)

########################
## Running the models ##
########################

## Stochastic
stoch_run <- ssa.adaptivetau(init.values = stoch_init.values, transitions = SIS_transitions, 
                       rateFunc = stoch_SIS_rateFunc, params = parameters, tf = 100)

stoch_run_df <- data.frame(stoch_run) # save stochastic model run as dataframe

## Deterministic
det_run <- ode(y = det_init.values, times = seq(0, 100, by = 1), func = det_SIS_rateFunc, 
               parms = parameters)

det_run_df <- as.data.frame(det_run) # change to data frame
det_run_df$time <- NULL # Delete time variable

############################################## Figures ######################################################

par(mfrow=c(1,2))

## Deterministic model
matplot(x = seq(0, 100, by = 1), y = det_run_df, type = "l",
        xlab = "Time", ylab = "Proportion susceptible or infected", main = "Deterministic SIS Model",
        lwd = 1, lty = 1, bty = "l", col = c("black","red"))

## Add legend
legend(60, 1.0, c("Susceptible", "Infected"), pch = 1, col = c("black","red"), bty = "n")

## Stochastic model
plot(x = stoch_run_df$time, y = stoch_run_df$I, type = "line", col = "red", ylim = c(0,10000),
     xlab = "Time", ylab = "Number susceptible or infected", main = "Stochastic SIS Model")
par(new=T)
plot(x = stoch_run_df$time, y = stoch_run_df$S, type = "line", ylim = c(0,1000), ylab = "") # susceptibles

## Add legend
legend(80, 10000, c("Susceptible", "Infected"), pch = 1, col = c("black", "red"), bty = "n")
