# MCMC estimating stochastic processes with a deterministic likelihood
# 3 unknowns: beta, gamma, and the time of infection for each individual
# Assuming that time of infection is deterministic
# 18/12/17

#######################
## Required packages ##
#######################
#install.packages("deSolve")
if (!require("deSolve")) install.packages("deSolve")
library("deSolve") #package for solving differential equations

###############
## Read data ##
###############
setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")
run_stoch <- read.csv("run_stoch.csv")

##################################
## Likelihood, prior, posterior ##
##################################

# Time
timestep <- run_stoch$time[2] - run_stoch$time[1]
end <- max(run_stoch$time)
times <- seq(0, end, by = timestep)

# Initial population: N-1 susceptible, 1 infectious, 0 recovered
init.values = c(
  S = run_stoch$S[1],
  I = run_stoch$I[1],
  R = run_stoch$R[1]
)
N = sum(init.values)

sir <- function(time, state, param) {
  
  # define model parameters in term of the natural parameters
  beta <- param[1,1] 
  gamma <- param[2,1]
  
  with(as.list(c(state, param)), {
    
    dS <- -beta * S * I 
    dI <- (beta * S * I) -(gamma * I)
    dR <-  gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Likelihood distribution
bg_likelihood <- function(param){
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  
  det_sir <- ode(y = init.values, times = times, func = sir, parms = param)
  det_sir <- as.data.frame(det_sir)
  
  total = array(0, dim = (c(nrow(run_stoch))))
  
  for (i in 1:nrow(run_stoch -1)){
    I = round(det_sir$I[i])
    S = (N - (I + run_stoch$R[i])) # Susceptibles for timestep i
    new_I = if (i == 1){
      I
    } else {
      (round(det_sir$I[i+1]) - I + run_stoch$R[i+1] - run_stoch$R[i]) # new I for timestep i+1
    }
  
    # print(c(S, new_I, det_sir$I[i], run_stoch$new_R[i], run_stoch$R[i], (S + I + run_stoch$R[i])))
    
    if (is.na(new_I) == F & (S<0 | new_I<0 | I<0)){ # (S<0 | new_I<0 | I<0)
      betalikelihood = -1000
      gammalikelihood = -1000
    } else {
      betalikelihood = dbinom(new_I, S, (1-(exp(-beta*I*timestep))), log = T)
      gammalikelihood = dbinom(run_stoch$new_R[i+1], I, (1-(exp(-gamma*timestep))), log = T)
    }
    
    # To deal with -Inf (i.e. log(0))
    if (is.na(new_I) == F & betalikelihood == -Inf){
      betalikelihood = -1000
    }
    if (is.na(new_I) == F & gammalikelihood == -Inf){
      gammalikelihood = -1000
    }
    total[i] = (betalikelihood + gammalikelihood)
  }
  
  ll = sum(total, na.rm = T)
  return(ll)
}

# Prior distribution
prior <- function(param){
  beta = as.numeric(param[1,1])
  gamma = as.numeric(param[2,1])
  
  betaprior = dunif(beta, min = 0, max = 100, log = T)
  gammaprior = dunif(gamma, min = 0, max = 100, log = T)
  
  return(betaprior + gammaprior)
}

# Posterior distribution
bg_posterior <- function(param){
  return (bg_likelihood(param) + prior(param))
}

# Proposal function for beta and gamma
proposalfunction <- function(param){
  param[1, 1] = rnorm(1, mean = param[1, 1], sd = 0.00004) # beta proposal rnorm(1, mean = 0.005, sd = 0.00001)
  param[2, 1] = rnorm(1, mean = param[2, 1], sd = 0.0004) # gamma proposal rnorm(1, mean = 0.08, sd = 0.00001)
  # print(c(param[1,1,1], param[2,1,1]))
  return(param)
}

##########
## MCMC ##
##########

run_metropolis_MCMC <- function(startvalue, iterations){
  divisor = divisor # the interval at which chain values are saved
  chain = array(dim = c((nrow(startvalue)+1), (iterations/divisor))) # Array for storing chain data
  temp_chain = array(dim = c((nrow(startvalue)+1), 2)) # Temporary array used for single iterations
  
  chain[1,1] = startvalue[1] # beta
  chain[2,1] = startvalue[2] # gamma
  chain[3,1] = 0             # likelihood

  
  temp_chain[,1] = chain[,1]
  
  for (i in 1:iterations){
    
    # Beta and gamma
    proposal = proposalfunction(temp_chain[1:2,]) # beta proposal and gamma proposal

    if (proposal[1,1] < 0.0 & proposal[2,1] < 0.0 ){
      temp_chain[1:2,2] = temp_chain[1:2,1]
      temp_chain[3,2] = temp_chain[3,1]
      
    } else if (proposal[1,1] < 0.0 & proposal[2,1] >= 0.0){
      proposal[1,1] = temp_chain[1,1]
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[1:2,])
      # print(c(bg_posterior(proposal), bg_posterior(temp_chain[,])))
      if (log(runif(1)) < probab){
        temp_chain[1:2,2] = proposal[1:2,1]
        temp_chain[3,2] = bg_likelihood(proposal)
      } else{
        temp_chain[1:2,2] = temp_chain[1:2,1]
        temp_chain[3,2] = bg_likelihood(temp_chain[1:2,])
      }
      
    } else if (proposal[1,1] >= 0.0 & proposal[2,1] < 0.0){
      proposal[2,1] = temp_chain[2,1]
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[1:2,])
      # print(c(bg_posterior(proposal), bg_posterior(temp_chain[,])))
      if (log(runif(1)) < probab){
        temp_chain[1:2,2] = proposal[1:2,1]
        temp_chain[3,2] = bg_likelihood(proposal)
      } else{
        temp_chain[1:2,2] = temp_chain[1:2,1]
        temp_chain[3,2] = bg_likelihood(temp_chain[1:2,])
      }
      
    } else{
      probab = bg_posterior(proposal) - bg_posterior(temp_chain[1:2,])
      # print(c(bg_posterior(proposal), bg_posterior(temp_chain[,])))
      if (log(runif(1)) < probab){
        temp_chain[1:2,2] = proposal[1:2,1]
        temp_chain[3,2] = bg_likelihood(proposal)
      } else{
        temp_chain[1:2,2] = temp_chain[1:2,1]
        temp_chain[3,2] = bg_likelihood(temp_chain[1:2,])
      }
    }
    
    # Save the iteration if it is divisible by the divisor without residuals
    if (i%%divisor == 0) {
      chain[,(i/divisor)] = temp_chain[,2]
      
    }    
    
    # Re-set temporary chain for next iteration
    temp_chain[,1] = temp_chain[,2]
    
    # Print every nth iteration, to know how far along run is
    if (i%%(iterations/20) == 0 | i == 1) {
      print(i)
      det_sir <- ode(y = init.values, times = times, func = sir, parms = temp_chain[1:2,])
      det_sir <- as.data.frame(det_sir)
      
      S = array(0, dim = (c(nrow(run_stoch))))
      new_I = array(0, dim = (c(nrow(run_stoch))))
      
      for (i in 1:nrow(run_stoch -1)){
        S[i] = (N - (round(det_sir$I[i]) + run_stoch$R[i])) # Susceptibles for timestep i
        new_I[i] = if (i == 1){
          round(det_sir$I[i])
        } else {
          (round(det_sir$I[i+1]) - round(det_sir$I[i]) + run_stoch$R[i+1] - run_stoch$R[i]) # new I for timestep i+1
        }
      }

      par(mfrow = c(2,1))
      
      plot(run_stoch$R, ylim = c(0, N), type = "l", col = "orange", xlab = "Timestep", ylab = "Number of individuals")
      lines(round(det_sir$I), type = "l", col = "red", xlab = " ", ylab = " ")
      lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")
      lines(round(det_sir$R), type = "l", col = "black", xlab = "", ylab = "")
      lines(S, type = "l", col = "darkolivegreen3", xlab = "", ylab = "")
      legend(100, 0.5*N, c("Deterministic recovered", "True recovered", "Deterministic infected", "True infected", "Susceptible"), pch = 1, col = c("black", "orange", "red", "grey", "darkolivegreen3"), bty = "n")
    
      plot(new_I, ylim = c(-10, N*0.25), type = "l", col = "red", xlab = "Timestep", ylab = "Number of individuals")
      # lines(run_stoch$new_R, type = "l", col = "orange", xlab = "", ylab = "")
      lines(run_stoch$new_I, type = "l", col = "grey", xlab = "", ylab = "")
      legend(100, 0.5*(N*0.25), c("Newly infected", "True newly infected"), pch = 1, col = c("red", "grey"), bty = "n")
      
    }
    
  }
  return(chain)
}

# Where to start the chain for beta, gamma, and I
startvalue <- array(dim = c(2))
startvalue[1] <-  beta # beta guess
startvalue[2] <-  gamma # gamma guess

# Number of runs
iterations =  5
divisor = 1 # how often runs are being saved

# Run the MCMC
set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)
rownames(chain) <- c("beta", "gamma", "loglik")

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.1*(iterations/divisor)
acceptance <- 1-mean(duplicated(chain[,-(1:burnIn)]))

################
## MCMC Plots ##
################

par(mfrow = c(2,2))

hist(chain[1,-(1:burnIn)],nclass=30, main="Posterior of beta")
abline(v = mean(chain[1,-(1:burnIn)]), col = "red")

hist(chain[2, -(1:burnIn)],nclass=30, main="Posterior of gamma")
abline(v = mean(chain[2,-(1:burnIn)]), col = "red")

plot(chain[1,], type = "l", main = "Chain values of beta")

plot(chain[2,], type = "l", main = "Chain values of gamma")

# Plot beta vs. gamma
par(mfrow = c(1,1))
library(RColorBrewer)
library(MASS)

plot(x = chain[2,], y = chain[1,], xlab = "Gamma", ylab = "Beta", pch = 20, cex = 0.8)

# k <- 11
# my.cols <- rev(brewer.pal(k, "RdYlBu"))
# z <- kde2d(chain[2,], chain[1,], n=50)
# filled.contour(z, nlevels=k, col=my.cols, xlab = "Gamma", ylab = "Beta")

#####################
## Likelihood plot ##
#####################

par(mfrow = c(1,2), mar=c(5,6,2,0.5))
plot(chain[3,], type = "l", main = "Chain values of log likelihood", xlab = "", ylab = "Log(likelihood)")
mtext("Iteration x100",side=1,line=2)
plot(chain[3,-(1:burnIn)], type = "l", main = "Zoomed in" , xlab = "", ylab = "Log(likelihood)")
mtext("(Iteration - burn-in) x100",side=1,line=2)

########################################################################################################################

########################################################################################################################

beta_gamma_loglik <- t(chain[,]) # transpose beta, gamma, and log likelihood data so that rows become columns etc.
colnames(beta_gamma_loglik) <- c("beta", "gamma", "loglik") # name the columns

setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")

# Beta, gamma, and log likelihood data
write.csv(data.frame(beta_gamma_loglik), file = "det_mcmc_beta_gamma_loglik_test.csv", row.names = FALSE)
