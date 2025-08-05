# Title: Stochastic simulation algorithm of the stage-structured mosquito pop. model
# Author: Joseph Baafi
# Date: 23/11/2021

# Clear workspace
rm(list = ls())

# Set working directory
setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/stochastic-simulations-R")

# Load packages
library(GillespieSSA)
library(ggplot2)

# The stage-structured mosquito population model is defined as:
# dE/dt = alpha*phi*A*(1-A\K) - (delta_E + mu_E)*E
# dL/dt = delta_E*E - (delta_L + mu_L)*L
# dP/dt = delta_L*L - (delta_P + mu_P)*P
# dA/dt = tau*delta_P*P - mu_A*A

#data = data.frame(out_tau.leap$data)
#write.csv(data, "time_df.csv", row.names = FALSE)

# data2 <- data.frame(out_tau.leap$data)
# write.csv(data2, "time_df2.csv", row.names = FALSE)
# t <- data$t

# Note::
# -----------------------------
# The parameters are functions of time and to be able to use them in the code, the
# time should be updated with the time used to run the stochastic model which is chosen 
# at random. I need  a code that is able to use a starting time, t=0 to run the model and update
# this time to run subsequent ones until the final time, ft=final is reached. 
# How do we do this in R? Ask Amy! 
# ------------------------------

# From Amy 
# 1. Differentiate temperature function with respect to time and put it into the
#    propensity function.
# 2. You can use a for-loop instead of ssa() to run the model. This way you will
#    be able to update the time for each time step (delta t) but that will be slow. 

# To-do (24/11/2021)
# 1. Study how to update the time from the stochastic simulation algorithm and use it
# in the next run
# 2. I can probably start with t=0 and use it to run the functions for tomorrow's meeting.
# If I start with zero, how do I update the subsequent t values?


# Define temperature function
temp <- function(t){
  temp = 8.923 - 8.764*sin(2 * pi * t/365) - 8.764*cos(2 * pi * t/365) + rnorm(n = 1, mean = 0, sd = 20.0) 
  return(temp)
}

plot(t, temp(t), "l")

# Functional oviposition rate
t <- data$t
phi <- function(t){
  func <- 6*exp(t/10)
  return(func)
}

#Define parameters 2 (source: Abdelrazec & Gumel, 2017)
parms2 <- c(
  alpha = 1,          # number of eggs laid
  phi = 10.7,        # oviposition
  K = 10^6,           # oviposition
  delta_E = 1/2.5,    # hatching rate of eggs into larvae
  delta_L = 1/7,      # development rate of larvae into pupae
  delta_P = 1/3,      # development rate of pupae into adult
  tau = 1/2,          # fraction of juviniles becoming adult females
  mu_E = 0.36,        # egg mortality rate
  mu_L = 0.34,        # larva mortality rate
  mu_P = 0.17,        # pupal mortality rate
  mu_A = 0.05         # adult mortality
)

tf <- 20                        # Final time
simName <- "Mosquito Population Model"   # Model name

# Define initial state vector
x0 <- c(E=10, L=1, P=1, A=1)

# Define state-change matrix
nu <- matrix(
  c(1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, -1),
  nrow = 4,
  byrow = TRUE)

# Define propensity functions
a <- c("alpha*phi*A*(1-(A/K))", "delta_E*E", "delta_L*L", "tau*delta_P*P",
       "mu_E*E", "mu_L*L", "mu_P*P", "mu_A*A")

# Run simulations with the Explict tau-leap method
out_etl <- ssa(x0 = x0, a = a, nu = nu, parms = parms2, tf = tf, method = ssa.etl(),
  simName = simName, verbose = FALSE, consoleInterval = 0)

ssa.plot(out_etl, show.title = TRUE, show.legend = TRUE)

# Plot
plot(out_etl$data[,1], out_etl$data[,2], col="red", cex=0.2, pch=19, 
     type = "l", main = "Mosquito Population", xlab = "Time (days)", ylab = "Population")
lines(out_etl$data[,1], out_etl$data[,3], col="green", cex=0.2, pch=19)
lines(out_etl$data[,1], out_etl$data[,4], col="blue", cex=0.2, pch=19)
lines(out_etl$data[,1], out_etl$data[,5], col="violet", cex=0.2, pch=19)
legend(x = "topleft",                                # Position
       legend = c("Egg", "Larva", "Pupa", "Adult"),  # Legend texts
       col = c("red", "green", "blue", "violet"),   # Line colors
       lwd = 2)

# Save the data into a data-frame
data = data.frame(out_etl$data)

# Plot with the ggplot function
ggplot(data=data, aes(x = t)) + 
  geom_line(aes(y = E), color = "red") +
  geom_line(aes(y = L), color = "green") +
  geom_line(aes(y = P), color = "blue") +
  geom_line(aes(y = A), color = "violet") +
  xlab('Time (days)') +
  ylab('Population')+
  theme_light() 









