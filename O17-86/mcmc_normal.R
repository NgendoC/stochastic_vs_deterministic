## MCMC for conjugate Normal distribution
# 19/10/2017

# Prior distiribution of theta has a mean of 0 and variance of 1
# Observation of data has an unknown mean and a variance of 1

# Packages needed
library("MCMCpack")

####################################################### 
## Use MCMCpack as a comparison for my manual method ##
#######################################################

# The data, known variance of data, prior mean of mu, prior variance of mu, number of MC draws
# don't know mean of data!
x <- seq(-2, 6, 0.01)
D <- dnorm(x, 2, 1) # A vector of observed data
mcmc_package <- MCnormalnormal(D, 1, 0, 1, mc = 1000)
plot(mcmc_package)
var(mcmc_package)
mean(mcmc_package)

######################
## My manual method ##
######################

# Use the same test data as for MCMCpack
plot(x, D, type='l', ylab="p(X)", xlab="x")

# Define the prior distribution (mean=0, variance=1)
mean_prior <- dnorm(x, 0, 1)

# Define the posterior distribution
mean_posterior <- (D * mean_prior)

plot(x, mean_posterior, type="l", ylim=c(0,0.6))
par(new=TRUE)
plot(x, mean_prior, type="l", col="red", ylim=c(0,0.6))

# Apply the Metropolis algorithm

# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
# proposalfunction <- function(param){
#   return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
# }
# run_metropolis_MCMC <- function(startvalue, iterations){
#   chain = array(dim = c(iterations+1,3))
#   chain[1,] = startvalue
#   for (i in 1:iterations){
#     proposal = proposalfunction(chain[i,])
#     
#     probab = exp(posterior(proposal) - posterior(chain[i,]))
#     if (runif(1) < probab){
#       chain[i+1,] = proposal
#     }else{
#       chain[i+1,] = chain[i,]
#     }
#   }
#   return(chain)
# }



