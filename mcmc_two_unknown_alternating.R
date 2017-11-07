## MCMC with two unknowns, alternating the change
## 31/10/2017

##################
## Known values ##
##################

true_mean <- 5
true_var <- 20
true_sd <- sqrt(true_var)
n_obs <- 31

# create values
D <- rnorm(n=n_obs,mean=true_mean,sd=true_sd)

##################################
## Likelihood, prior, posterior ##
##################################

mean_likelihood <- function(param){
  mean = param
  singlelikelihoods = dnorm(D, mean = mean, sd = sqrt(1), log = T)
  sumll = sum(singlelikelihoods)
  return(sumll) 
}

var_likelihood <- function(param){
  var = param
  singlelikelihoods = dnorm(D, mean = 1, sd = sqrt(var), log = T)
  sumll = sum(singlelikelihoods)
  return(sumll) 
}

# Prior distribution
mean_prior <- function(param){
  mean = param
  meanprior = dnorm(mean, mean = 0, sd = 1, log = T)
  return(meanprior)
}

var_prior <- function(param){
  var = param
  varprior = dnorm(var, mean = 1, sd = 1, log = T)
  return(varprior)
}

# Posterior distribution
mean_posterior <- function(param){
  return (mean_likelihood(param) + mean_prior(param))
}

var_posterior <- function(param){
  return (var_likelihood(param) + var_prior(param))
}

##########
## MCMC ##
##########

proposalfunction <- function(param){
  return(rnorm(1,mean = param, sd= 1))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  
  for (i in 1:iterations){
    # Chain for the mean
    mean_proposal = proposalfunction(chain[i,1])
    mean_probab = exp(mean_posterior(mean_proposal) - mean_posterior(chain[i,1]))
    if (runif(1) < mean_probab){
      chain[i+1,1] = mean_proposal
    }else{
      chain[i+1,1] = chain[i,1]
    }
    
    # Chain for the variance
    var_proposal = proposalfunction(chain[i,2])
    var_probab = exp(var_posterior(var_proposal) - var_posterior(chain[i,2]))
    if (runif(1) < var_probab){
    chain[i+1,2] = var_proposal
    }else{
    chain[i+1,2] = chain[i,2]
    }
  }
  
  return(chain)
}

# Where to start the chain
# Takes a random number from the (Normal) prior distribution
startvalue <- c(1,1) #rnorm(2, prior_mean, prior_sd)

# Number of runs
iterations = 10000
set.seed(4)
chain <- run_metropolis_MCMC(startvalue, iterations)

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.1*iterations
acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

###########
## Plots ##
###########

par(mfrow = c(2,2))

hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of mean")
abline(v = mean(chain[-(1:burnIn),1]), col = "red")

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of var")
abline(v = mean(chain[-(1:burnIn),2]), col = "red")

plot(chain[-(1:burnIn),1], type = "l", main = "Chain values of mean")

plot(chain[-(1:burnIn),2], type = "l", main = "Chain values of var")

###################################################################################################################

# likelihood <- function(param){
#   mean = param[1]
#   var = param[2]
#   
#   singlelikelihoods = dnorm(D, mean = mean, sd = sqrt(var), log = T)
#   sumll = sum(singlelikelihoods)
#   return(sumll) 
# }

# prior <- function(param){
#   mean = param[1]
#   var = param[2]
#   meanprior = dnorm(mean, mean = 0, sd = 1, log = T)
#   varprior = dnorm(var, mean = 0, sd = 1, log = T)
#   return(meanprior+varprior)
# }

# posterior <- function(param){
#   return (likelihood(param) + prior(param))
# }


# Helpful resource for MCMC in R:
# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
