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

# Make prior distribution of theta have a mean = 0 and variance = 1
set.seed(1)
xseq <- seq(-4, 4, .01)
prior <- dnorm(xseq, 0, 1)
plot(xseq, prior, xlim=c(-4, 4), type='l', ylab="p(X)", xlab="x")

# Make observation mean = 2 and variance = 1 - in real life I wouldn't know the mean but I need to make data
D <-  dnorm(xseq, 2, 1)
plot(xseq, D, xlim=c(-4, 4), type='l', ylab="p(X)", xlab="x")

# Plot both the prior and the observation
plot(xseq, prior, xlim=c(-4, 4), type='l', ylab="p(X)", xlab="x")
par(new=T)
plot(xseq, D, xlim=c(-4, 4), type='l', ylab="", xlab="", col = "red")

