## MCMC for conjugate Normal distribution
# 19/10/2017

# Prior distiribution of theta has a mean of 0 and variance of 1
# Observation of data has an unknown mean and a variance of 1

####################################################### 
## Use MCMCpack as a comparison for my manual method ##
#######################################################

# Test my manual method against this
library("MCMCpack")
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





