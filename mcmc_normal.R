## MCMC for conjugate Normal distribution
# 19/10/2017

# Prior distiribution of theta has a mean of 0 and variance of 1
# Observation of data has an unknown mean and a variance of 1

# Packages needed
library("MCMCpack")

################# 
## Create data ##
#################

obs_var <- 1
obs_sd <- obs_var^(1/2)
obs_mean <- 2

prior_var <- 1
prior_sd <- prior_var^(1/2)
prior_mean <- 0

# create independent x-values 
x <- seq(-2, 6, 0.01)
# create dependent values according to Normal distribution
D <- dnorm(x, obs_mean, obs_var)

plot(x,D, main="Test Data")

####################################################### 
## Use MCMCpack as a comparison for my manual method ##
#######################################################

# The data, known variance of data, prior mean of mu, prior variance of mu, number of MC draws
# don't know mean of data!
mcmc_package <- MCnormalnormal(D, obs_var, prior_mean, prior_var, mc = 1000)
plot(mcmc_package)
var(mcmc_package)
mean(mcmc_package)

######################
## My manual method ##
######################

# Use the same test data as for MCMCpack
# This may be wrong
likelihood <- dnorm(x, obs_mean, obs_sd, log= T)

# Define the prior distribution (mean=0, variance=1)
mean_prior <- dnorm(x, prior_mean, prior_sd, log=T)

# Define the posterior distribution
# Note: is on a logarithmic scale
mean_posterior <- likelihood + mean_prior

plot(x, mean_posterior, type="l", col = "red", ylim = c(-30, 0), ylab = "log(mean)")
par(new=TRUE)
plot(x, mean_prior, type="l", col="black", ylim = c(-30, 0), ylab = " ")

# Apply the Metropolis algorithm
# Start at random parameter value
# Choose new parameter close to old value based on prob. density ratio
# Jump to new probability if P(old)/P(new) >1


# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
# This example works in logs!
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

# metropolis=function(n=1000,eps=5) 
# {
#   vec=vector("numeric", n)
#   x=0
#   vec[1]=x
#   for (i in 2:n) {
#     innov=runif(1,-eps,eps)
#     can=x+innov
#     aprob=min(1,dnorm(can)/dnorm(x))
#     u=runif(1)
#     if (u < aprob) 
#       x=can
#     vec[i]=x
#   }
#   vec
# }



