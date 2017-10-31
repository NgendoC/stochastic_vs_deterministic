## MCMC with two unknowns
## 31/10/2017

##################
## Known values ##
##################

trueMean <- 5
trueSd <- 10
sampleSize <- 31

# create values according to a + N(0,sd)
y <- rnorm(n=sampleSize,mean=trueMean,sd=trueSd)

##################################
## Likelihood, prior, posterior ##
##################################

likelihood <- function(param){
  mean = param[1]
  sd = param[2]
  
  singlelikelihoods = dnorm(y, mean = mean, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

# Prior distribution
prior <- function(param){
  mean = param[1]
  sd = param[2]
  meanprior = dnorm(mean, mean = 0, sd = 1, log = T)
  sdprior = dnorm(sd, mean = 0, sd = 1, log = T)
  return(meanprior+sdprior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

##########
## MCMC ##
##########

proposalfunction <- function(param){
  return(rnorm(2,mean = param, sd= c(0.1,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

# Number of runs
iterations = 10000
set.seed(4)
chain <- metropolis_MCMC(startvalue, iterations)
# 
# # The beginning of the chain is biased towards the starting point, so take them out
# # normally burnin is 10%-50% of the runs
burnIn = 0.1*iterations
acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

###########
## Plots ##
###########

par(mfrow = c(2,2))

hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of mean")
abline(v = mean(chain[-(1:burnIn),1]), col = "red")

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of sd")
abline(v = mean(chain[-(1:burnIn),2]), col = "red")

plot(chain[-(1:burnIn),1], type = "l", main = "Chain values of a")

plot(chain[-(1:burnIn),2], type = "l", main = "Chain values of sd")

###################################################################################################################
## Previous attempt
# 
# ################ 
# ## Known data ##
# ################
# 
# # Data observation
# #D <- 2 # 1 observation
# D <- rnorm(30, mean = 2, sd = 1)  # many observations, random numbers from a normal distribution
# 
# # Number of observations
# num_obs <- length(D)
# 
# # Mean and variance for the prior distribution 
# prior_var <- 1
# prior_sd <- sqrt(prior_var)
# prior_mean <- 0
# 
# # Proposal variance
# proposal_var <- 1
# proposal_sd <- sqrt(proposal_var)
# 
# #################
# ## Manual MCMC ##
# #################
# 
# # Likelihood function for the observation
# likelihood <- function(param){
#   mean = param[1]
#   var = param[2]
#   
#   obs_likelihood = dnorm(D, mean, sqrt(var), log = T)
#   sum_likelihood = sum(obs_likelihood)
#   return(sum_likelihood)
# } 
# 
# # Define the prior distribution
# prior <- function(param){
#   mean = param[1]
#   var = param[2]
#   
#   mean_prior = dnorm(mean, mean = 0, sd = 1, log = T)
#   var_prior = dnorm(var, mean = 0, sd = 1, log = T)
#   return(mean_prior + var_prior)
# }
# 
# # Define the posterior distribution
# # Note: is on a logarithmic scale
# posterior <- function(param){
#   return (likelihood(param) + prior(param))
# }
# 
# # Apply the Metropolis algorithm
# 
# # Probability density distribution for proposing new values
# # Proposal function has a Normal distribution 
# # The distribution is centered around the current value in the Markov chain
# proposalfunction <- function(param){
#   return(rnorm(2, param, sd = c(1,1)))
# }
# 
# # Start at random parameter value
# # Choose new parameter close to old value based on prob. density ratio
# 
# metropolis_MCMC <- function(startvalue, iterations){
#   chain = array(dim = c(iterations+1,2))
#   
#   # Give a starting point to the chain
#   chain[1,] = startvalue
#   
#   # For the other iterations, use the proposal function to suggest the next value
#   for (i in 1:iterations){
#     proposal = proposalfunction(chain[i,])
#     
#     # P(new)/P(old) NOTE: LOGARITHMIC
#     prob = exp(posterior(proposal) - posterior(chain[i,]))
#     
#     # Jump to new probability if P(new)/P(old) >1
#     # if r < than P, accept, else reject
#     if (runif(1) < prob){ 
#       chain[i+1,] = proposal
#     }else{
#       chain[i+1,] = chain[i,]
#     }
#   }
#   return(chain)
# }
# 
# # Where to start the chain
# # Takes a random number from the (Normal) prior distribution
# startvalue <- c(1,1) #rnorm(2, prior_mean, prior_sd)
# 
# # Number of runs
# iterations = 10000
# set.seed(4)
# chain <- metropolis_MCMC(startvalue, iterations)
# 
# # The beginning of the chain is biased towards the starting point, so take them out
# # normally burnin is 10%-50% of the runs
# burnIn = 0.1*iterations
# acceptance <- 1-mean(duplicated(chain[-(1:burnIn)]))
# 
# ###########
# ## Plots ##
# ###########
# 
# # Plot the MCMC
# par(mfrow = c(1,2))
# hist(chain[-(1:burnIn)],nclass=30, main="Posterior mean", freq=FALSE, xlim=c(-2,4), ylim=c(0,1.0), col = "grey")
# abline(v = mean(chain[-(1:burnIn)]))
# plot(chain[-(1:burnIn)], type = "l", main = "Chain values of mean")


# Helpful resource for MCMC in R:
# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/