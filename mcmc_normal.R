## MCMC for conjugate Normal distribution
# 19/10/2017

# Prior distiribution of theta has a mean of 0 and variance of 1
# Observation of data has an unknown mean and a variance of 1

################ 
## Known data ##
################

# Data observation and variance
# Note: mean is not known!
#D <- 2 # 1 observation
D <- rnorm(4, mean = 2, obs_sd)  # many observations, random numbers from a normal distribution

# Number of observations
num_obs <- length(D)
mean_D <- mean(D)

# Mean and variance for the observed data
obs_var <- 1
obs_sd <- sqrt(obs_var)
num_obs_sd <- sqrt((num_obs/obs_var)^-1)

# Mean and variance for the prior distribution of the unknown mean
prior_var <- 1
prior_sd <- sqrt(prior_var)
prior_mean <- 0

# Calculation for true posterior mean distribution
# mean
true_mean <- (((obs_var/num_obs)*prior_mean) + (prior_var*mean_D)) / ((obs_var/num_obs) + prior_var)
true_var <- ((num_obs/obs_var) + (1/prior_var))^-1
true_sd <- sqrt(true_var)

# Proposal variance
proposal_var <- 1
proposal_sd <- sqrt(proposal_var)

####################################################### 
## Use MCMCpack as a comparison for my manual method ##
#######################################################
# Packages needed
#library("MCMCpack")

# # Create independent x-values 
# x <- seq(-2, 6, 0.01)
# # Create dependent values according to Normal distribution
# D <- dnorm(x, obs_mean, obs_var)
# 
# plot(x,D, main="Test Data", type="l")

# The data, known variance of data, prior mean of mu, prior variance of mu, number of MC draws
# don't know mean of data!
# set.seed(4)
# mcmc_package <- MCnormalnormal(D, obs_var, prior_mean, prior_var, mc = 1000)
# plot(mcmc_package)
# var(mcmc_package)
# mean(mcmc_package)

######################
## My manual method ##
######################

# Likelihood function for the observation
likelihood <- function(param){
  obs_likelihood = dnorm(D, param, num_obs_sd, log= T)
  return(obs_likelihood)
} 

# Define the prior distribution (mean=0, variance=1)
mean_prior <- function(param){
  prior = dnorm(param, prior_mean, prior_sd, log = T)
  return(prior)
}

# Define the posterior distribution
# Note: is on a logarithmic scale
mean_posterior <- function(param){
  return (likelihood(param) + mean_prior(param))
}

# Apply the Metropolis algorithm

# Probability density distribution for proposing new values
# Proposal function has a Normal distribution 
# The distribution is centered around the current value in the Markov chain
proposalfunction <- function(param){
  return(rnorm(1, param, proposal_sd))
}

# Start at random parameter value
# Choose new parameter close to old value based on prob. density ratio

metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1))
  
  # Give a starting point to the chain
  chain[1] = startvalue
  
  # For the other iterations, use the proposal function to suggest the next value
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i])
    
    # P(new)/P(old) NOTE: LOGARITHMIC
    prob = exp(mean_posterior(proposal) - mean_posterior(chain[i]))
    
    # Jump to new probability if P(new)/P(old) >1
    # if r < than P, accept, else reject
    if (runif(1) < prob){ 
      chain[i+1] = proposal
    }else{
      chain[i+1] = chain[i]
    }
  }
  return(chain)
}

# Where to start the chain
# Takes a random number from the (Normal) prior distribution
startvalue <- rnorm(1, prior_mean, prior_sd)

# Number of runs
iterations = 10000
set.seed(4)
chain <- metropolis_MCMC(startvalue, iterations)

# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
burnIn = 0.1*iterations
acceptance <- 1-mean(duplicated(chain[-(1:burnIn)]))

#####################################
## Outcomes from the manual method ##
#####################################

# Plot the MCMC
par(mfrow = c(1,2))
hist(chain[-(1:burnIn)],nclass=30, main="Posterior mean", freq=FALSE, xlim=c(-2,4), ylim=c(0,1.0), col = "grey")
abline(v = mean(chain[-(1:burnIn)]))
par(new=T)
xdata = seq(-2, 4, 0.1)
plot(xdata, dnorm(xdata, true_mean, true_sd), type = "l", xlab = " ", ylab = " ",
     xlim=c(-2,4), ylim=c(0,1.0), col="red", lty=2, lwd=2)
plot(chain[-(1:burnIn)], type = "l", main = "Chain values of mean")


# Helpful resource for MCMC in R:
# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
