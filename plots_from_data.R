# Making plots from saved data

setwd("C:/Users/Janetta Skarp/OneDrive - Imperial College London/MRes_BMR/Project_1/Work_folder/Data")


re_bootstrap <- read.csv("re_betagamma_bootstrap_27.01.18.csv")
re_point <- read.csv("re_betagamma_point_27.01.18.csv")
  
stoch_mcmc <- read.csv("stoch_mcmc_betagamma_25.01.18.csv")

det_mcmc <- read.csv("det_mcmc_betagamma_25.01.18.csv")

burnIn = 250

########################
## With constant xlim ##
########################
# Histogram
par(mfrow = c(2,3))

# RE beta
hist(as.numeric(re_bootstrap[1,2:ncol(re_bootstrap)]),nclass=30, main="RE Beta", xlab="Beta", xlim= c(0.002, 0.007))
abline(v = re_point[1,2], col = "red")

# Det MCMC beta
hist(as.numeric(det_mcmc[1,burnIn:ncol(det_mcmc)]),nclass=30, main="Det MCMC Beta", xlab="Beta", xlim= c(0.002, 0.007))
abline(v = mean(as.numeric(det_mcmc[1,burnIn:ncol(det_mcmc)])), col = "red")

# Stoch MCMC beta
hist(as.numeric(stoch_mcmc[1,burnIn:ncol(stoch_mcmc)]),nclass=30, main="Stoch MCMC Beta", xlab="Beta", xlim= c(0.002, 0.007))
abline(v = mean(as.numeric(stoch_mcmc[1,burnIn:ncol(stoch_mcmc)])), col = "red")

# RE gamma
hist(as.numeric(re_bootstrap[2,2:ncol(re_bootstrap)]),nclass=30, main="RE Gamma", xlab="Gamma", xlim=c(0.05,0.11))
abline(v = re_point[2,2], col = "red")

# Det MCMC gamma
hist(as.numeric(det_mcmc[2, burnIn:ncol(det_mcmc)]),nclass=30, main="Det MCMC Gamma", xlab="Gamma", xlim=c(0.05,0.11))
abline(v = mean(as.numeric(det_mcmc[2,burnIn:ncol(det_mcmc)])), col = "red")

# Stoch MCMC gamma
hist(as.numeric(stoch_mcmc[2, burnIn:ncol(stoch_mcmc)]),nclass=30, main="Stoch MCMC Gamma", xlab="Gamma", xlim=c(0.05,0.11))
abline(v = mean(as.numeric(stoch_mcmc[2,burnIn:ncol(stoch_mcmc)])), col = "red")

##################
## Without xlim ##
##################
# Histogram
par(mfrow = c(2,3))

# RE beta
hist(as.numeric(re_bootstrap[1,2:ncol(re_bootstrap)]),nclass=30, main="RE Beta", xlab="Beta")
abline(v = re_point[1,2], col = "red")

# Det MCMC beta
hist(as.numeric(det_mcmc[1,250:ncol(det_mcmc)]),nclass=30, main="Det MCMC Beta", xlab="Beta")
abline(v = mean(as.numeric(det_mcmc[1,250:ncol(det_mcmc)])), col = "red")

# Stoch MCMC beta
hist(as.numeric(stoch_mcmc[1,250:ncol(stoch_mcmc)]),nclass=30, main="Stoch MCMC Beta", xlab="Beta")
abline(v = mean(as.numeric(stoch_mcmc[1,250:ncol(stoch_mcmc)])), col = "red")

# RE gamma
hist(as.numeric(re_bootstrap[2,2:ncol(re_bootstrap)]),nclass=30, main="RE Gamma", xlab="Gamma")
abline(v = re_point[2,2], col = "red")

# Det MCMC gamma
hist(as.numeric(det_mcmc[2, 250:ncol(det_mcmc)]),nclass=30, main="Det MCMC Gamma", xlab="Gamma")
abline(v = mean(as.numeric(det_mcmc[2,250:ncol(det_mcmc)])), col = "red")

# Stoch MCMC gamma
hist(as.numeric(stoch_mcmc[2, 250:ncol(stoch_mcmc)]),nclass=30, main="Stoch MCMC Gamma", xlab="Gamma")
abline(v = mean(as.numeric(stoch_mcmc[2,250:ncol(stoch_mcmc)])), col = "red")

###################################
## Plot code from residual error ##
###################################

re_data1 <- read.csv("re_betagamma_badstart_01.02.18.csv") # Starting parameter guess not based on point estimate
re_data2 <- read.csv("re_betagamma_goodstart_01.02.18.csv") # Starting parameter guess based on point estimate

# print(c(re_data$beta[1], re_data$gamma[1]))

# Histogram
par(mfrow = c(1,3))

# Beta
hist(re_data1$beta[2:nrow(re_data1)],nclass=30, col = rgb(0.1,0.1,0.1,0.5), main="Beta", xlab="Beta value")
abline(v = re_data1$beta[1], col = "red")
hist(re_data2$beta[2:nrow(re_data2)],nclass=30, col=rgb(0.8,0.8,0.8,0.5), add = T)
abline(v = re_data2$beta[1], col = "red")
box()

# Gamma
hist(re_data1$gamma[2:nrow(re_data1)],nclass=30, col = rgb(0.1,0.1,0.1,0.5), main="Gamma", xlab="Gamma value")
abline(v = re_data1$gamma[1], col = "red")
hist(re_data2$gamma[2:nrow(re_data2)],nclass=30, col=rgb(0.8,0.8,0.8,0.5), add = T)
abline(v = re_data2$gamma[1], col = "red")
box()

# Residual error
hist(re_data1$RE[2:nrow(re_data1)],nclass=30, col=rgb(0.1,0.1,0.1,0.5), main="RE", xlab="Residual Error")
abline(v = re_data1$RE[1], col = "red")
hist(re_data2$RE[2:nrow(re_data2)],nclass=30, col=rgb(0.8,0.8,0.8,0.5), add = T)
abline(v = re_data2$RE[1], col = "red")
box()

# Beta vs. Gamma
par(mfrow = c(1,1))
plot(x = re_data1$gamma, y = re_data1$beta, col = rgb(1,0,0,0.5), xlab = "Gamma", ylab = "Beta", pch = 4, cex = 1)
  points(x = re_data2$gamma, y = re_data2$beta, col = rgb(0,0,0,0.5), xlab = "Gamma", ylab = "Beta", pch = 1, cex = 1)

# # Lines
# run_det <- as.data.frame(ode(y = init.values, times = times, func = sir, parms = sse_fit$par))
# 
# par(mfrow = c(1,1))
# plot(run_stoch$R, ylim = c(0, N), type = "l", col = "orange", xlab = "Timestep", ylab = "Number of individuals")
# lines(run_det$I, type = "l", col = "red", xlab = " ", ylab = " ")
# lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")
# lines(run_det$R, type = "l", col = "black", xlab = "", ylab = "")
# legend(100, 0.5*N, c("Deterministic recovered", "True recovered", "Deterministic infected", "True infected"), pch = 1, col = c("black", "orange", "red", "grey"), bty = "n")

# Making a heatmap for beta vs. gamma residual errors
if (!require("plotly")) install.packages("plotly")
library("plotly") #package for solving differential equations

# heatmap <- read.csv("re_heatmap_test.csv")
heatmap <- read.csv("re_heatmap_small_range.csv")
matrix_heatmap <- xtabs(RE~beta+gamma, data=heatmap)


beta <- seq(min(heatmap$beta), max(heatmap$beta), by = ((max(heatmap$beta)-min(heatmap$beta))/nrow(matrix_heatmap)))
gamma <- seq(min(heatmap$gamma), max(heatmap$gamma), by = ((max(heatmap$gamma)-min(heatmap$gamma))/nrow(matrix_heatmap)))

tick_beta <- list(
  autotick = FALSE,
  ticks = "outside",
  tick0 = min(heatmap$beta),
  dtick = ((max(heatmap$beta)-min(heatmap$beta))/nrow(matrix_heatmap)),
  ticklen = 5,
  tickwidth = 2,
  tickcolor = toRGB("blue")
)

tick_gamma <- list(
  autotick = FALSE,
  ticks = "outside",
  tick0 = min(heatmap$gamma),
  dtick = ((max(heatmap$gamma)-min(heatmap$gamma))/nrow(matrix_heatmap)),
  tickangle = 45,
  ticklen = 5,
  tickwidth = 2,
  tickcolor = toRGB("blue")
)

m <- list(
  l = 100,
  r = 5,
  b = 80,
  t = 5,
  pad = 4
)

# vals <- unique(scales::rescale(c(matrix_heatmap)))
# o <- order(10*vals, decreasing = FALSE)
# cols <- scales::col_numeric("Blues", domain = NULL)(-10*vals)
# colz <- setNames(data.frame(vals[o], cols[o]), NULL)
# p <- plot_ly(z = volcano, colorscale = colz, type = "heatmap")
# col_num = 80
# grays <- array(dim = c(col_num))
# for (i in 1:col_num){
#   grays[i] = paste("gray",i+(100-col_num), sep = "")
# }

plot_heatmap <- plot_ly(z = matrix_heatmap, x = ~gamma, y = ~beta, colors = colorRamp(c("yellow","darkorange2","orangered","red","maroon","magenta4","blue", "navy", "midnightblue")), type = "heatmap") %>%
  layout(xaxis = tick_gamma, yaxis = tick_beta, margin = m)
plot_heatmap

# rev(c("white", "yellow", "gold", "goldenrod1","orange", "darkorange","darkorange2", "orangered","red","firebrick3","firebrick","firebrick4", "deeppink4", "darkmagenta", "darkorchid4", "darkslateblue", "dodgerblue4", "dodgerblue3", "deepskyblue3", "turquoise3", "turquoise", "palegreen2", "palegreen3", "palegreen4"))),
# colorRamp(c("yellow", "darkorange2","orangered","red", "maroon","magenta4","blue", "navy", "midnightblue")),

####################################
## Plot code from stochastic MCMC ##
####################################

# plot(run_stoch$time, run_stoch$R, ylim = c(0,N), type = "l", col = "orange", xlab = "time (days)", ylab = "Number infectious/recovered")
# par(new=T)
# plot(run_stoch$time, run_stoch$guess_I, ylim = c(0,N), type = "l", col = "red", xlab = " ", ylab = " ")
# par(new=T)
# plot(x = run_stoch$time, y = run_stoch$I, type = "l", col = "black", ylim = c(0,N), xlab = " ", ylab = " ")
# legend(60, 0.8*N, c("Recovered", "Guessed infected", "True infected"), pch = 1, col = c("orange", "red", "black"), bty = "n")
#
# plot(run_stoch$guess_I, ylim = c(0, N), type = "l", col = "red", xlab = "Timestep", ylab = "Number of individuals infected")
# lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")
# lines(temp_chain[,1,2], type = "l", lty = 2, col = "black", xlab = " ", ylab = " ")
# legend(130, 1.0*N, c("True infected", "Guessed infected", "MCMC"), pch = 1, col = c("grey", "red", "black"), bty = "n")
#
# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
# burnIn = 0.1*(iterations/divisor)
# acceptance <- 1-mean(duplicated(chain[,-(1:burnIn),]))
# inf_acceptance <- 1-mean(duplicated(chain[,-(1:burnIn),2]))
#
# #Histogram
# par(mfrow = c(2,2))
# 
# hist(chain[1,-(1:burnIn),1],nclass=30, main="Posterior of beta")
# abline(v = mean(chain[1,-(1:burnIn),1]), col = "red")
# 
# hist(chain[2, -(1:burnIn),1],nclass=30, main="Posterior of gamma")
# abline(v = mean(chain[2,-(1:burnIn),1]), col = "red")
# 
# plot(chain[1, -(1:burnIn),1], type = "l", main = "Chain values of beta")
# 
# plot(chain[2, -(1:burnIn),1], type = "l", main = "Chain values of gamma")
# 
# # Plot beta vs. gamma
# par(mfrow = c(1,1))
# library(RColorBrewer)
# library(MASS)
# 
# plot(x = chain[2,,1], y = chain[1,,1], xlab = "Gamma", ylab = "Beta", pch = 20, cex = 0.8)
# 
# k <- 11
# my.cols <- rev(brewer.pal(k, "RdYlBu"))
# z <- kde2d(chain[2,,1], chain[1,,1], n=50)
# filled.contour(z, nlevels=k, col=my.cols, xlab = "Gamma", ylab = "Beta")
#  
# par(mfrow = c(1,1))
#   
# plot(run_stoch$guess_I, ylim = c(0, N), type = "l", col = "red", xlab = "Timestep", ylab = "Number of individuals infected")
#   lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")
#   lines(chain[,ncol(chain),2], type = "l", lty = 2, col = "black", xlab = " ", ylab = " ")
#   legend(130, 1.0*N, c("True infected", "Guessed infected", "MCMC"), pch = 1, col = c("grey", "red", "black"), bty = "n")
# 
# plot(run_stoch$guess_I, ylim = c(0, N), type = "l", col = "red", xlab = "Timestep", ylab = "Number of individuals infected")
#   lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")  
#   for (i in 1:ncol(chain)){
#     lines(chain[,i,2], type = "l", lty = 2, col = "black", xlab = " ", ylab = " ")
#   }
#   legend(130, 1.0*N, c("True infected", "Guessed infected", "MCMC"), pch = 1, col = c("grey", "red", "black"), bty = "n")  

#######################################
## Plot code from deterministic MCMC ##
#######################################

# det_sir <- ode(y = init.values, times = times, func = sir, parms = temp_chain[1:2,])
# det_sir <- as.data.frame(det_sir)
# 
# S = array(0, dim = (c(nrow(run_stoch))))
# new_I = array(0, dim = (c(nrow(run_stoch))))
# 
# for (i in 1:nrow(run_stoch -1)){
#   S[i] = (N - (round(det_sir$I[i]) + run_stoch$R[i])) # Susceptibles for timestep i
#   new_I[i] = if (i == 1){
#     round(det_sir$I[i])
#   } else {
#     (round(det_sir$I[i+1]) - round(det_sir$I[i]) + run_stoch$R[i+1] - run_stoch$R[i]) # new I for timestep i+1
#   }
# }
#
# par(mfrow = c(2,1))
# 
# plot(run_stoch$R, ylim = c(0, N), type = "l", col = "orange", xlab = "Timestep", ylab = "Number of individuals")
# lines(round(det_sir$I), type = "l", col = "red", xlab = " ", ylab = " ")
# lines(run_stoch$I, type = "l", col = "grey", xlab = " ", ylab = " ")
# lines(round(det_sir$R), type = "l", col = "black", xlab = "", ylab = "")
# lines(S, type = "l", col = "darkolivegreen3", xlab = "", ylab = "")
# legend(100, 0.5*N, c("Deterministic recovered", "True recovered", "Deterministic infected", "True infected", "Susceptible"), pch = 1, col = c("black", "orange", "red", "grey", "darkolivegreen3"), bty = "n")
# 
# plot(new_I, ylim = c(-10, N*0.25), type = "l", col = "red", xlab = "Timestep", ylab = "Number of individuals")
# # lines(run_stoch$new_R, type = "l", col = "orange", xlab = "", ylab = "")
# lines(run_stoch$new_I, type = "l", col = "grey", xlab = "", ylab = "")
# legend(100, 0.5*(N*0.25), c("Newly infected", "True newly infected"), pch = 1, col = c("red", "grey"), bty = "n")
# The beginning of the chain is biased towards the starting point, so take them out
# normally burnin is 10%-50% of the runs
# burnIn = 0.1*(iterations/divisor)
# acceptance <- 1-mean(duplicated(chain[,-(1:burnIn)]))
#
## MCMC Plots
# par(mfrow = c(2,2))
# 
# hist(chain[1,-(1:burnIn)],nclass=30, main="Posterior of beta")
# abline(v = mean(chain[1,-(1:burnIn)]), col = "red")
# 
# hist(chain[2, -(1:burnIn)],nclass=30, main="Posterior of gamma")
# abline(v = mean(chain[2,-(1:burnIn)]), col = "red")
# 
# plot(chain[1,], type = "l", main = "Chain values of beta")
# 
# plot(chain[2,], type = "l", main = "Chain values of gamma")
# 
# # Plot beta vs. gamma
# par(mfrow = c(1,1))
# library(RColorBrewer)
# library(MASS)
# 
# plot(x = chain[2,], y = chain[1,], xlab = "Gamma", ylab = "Beta", pch = 20, cex = 0.8)
#
# k <- 11
# my.cols <- rev(brewer.pal(k, "RdYlBu"))
# z <- kde2d(chain[2,], chain[1,], n=50)
# filled.contour(z, nlevels=k, col=my.cols, xlab = "Gamma", ylab = "Beta")
#
## Likelihood plot
# par(mfrow = c(1,2), mar=c(5,6,2,0.5))
# plot(chain[3,], type = "l", main = "Chain values of log likelihood", xlab = "", ylab = "Log(likelihood)")
# mtext("Iteration x100",side=1,line=2)
# plot(chain[3,-(1:burnIn)], type = "l", main = "Zoomed in" , xlab = "", ylab = "Log(likelihood)")
# mtext("(Iteration - burn-in) x100",side=1,line=2)
