# Title: Stochastic simulation algorithm of the stage-structured mosquito pop. model
# Author: Joseph Baafi
# Date: 18/11/2021

# Clear workspace
rm(list = ls())

# Set working director
#setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/stochastic-simulations-R")

# Load packages
library(GillespieSSA)
library(ssar)
library(ggplot2)

# The stage-structured mosquito population model is defined as:

# dE/dt = alpha*phi*A*(1-A\K) - (delta_E + mu_E)*E
# dL/dt = delta_E*E - (delta_L + mu_L + chi_L*L)*L
# dP/dt = delta_L*L - (delta_P + mu_P)*P
# dA/dt = tau*delta_P*P - mu_A*A

# #Define parameters (source: https://www.pascomosquito.org/resources/mosquito-biology/, https://megacatch.com/mosquito-lifecycle-faqs/)
# parms <- c(
#   alpha = 50,      # average number of eggs laid per oviposition (source: CDC)
#   phi = 1/4,
#   K = 10^6,          #0.051 # oviposition rate per adult mosquito (source: https://www.pascomosquito.org/resources/mosquito-biology/, https://megacatch.com/mosquito-lifecycle-faqs/)
#   delta_E = 1/2,    # hatching rate of eggs into larvae
#   delta_L = 1/4,    # development rate of larvae into pupae
#   delta_P = 1/2,    # development rate of pupae into adult
#   tau = 1/2,          # fraction of juviniles becoming adult females
#   mu_E = 0.36,      # egg mortality rate (SOurce: Abdelrazec & Gumel, 2017)
#   mu_L = 0.34,      # larva mortality rate
#   mu_P = 0.17,      # pupal mortality rate
#   mu_A = 0.05       # adult mortality
#   )

#Define parameters 2 (source: Abdelrazec & Gumel, 2017)
parms2 <- c(
  alpha = 1,      # number of eggs laid
  phi = 10.7,
   K = 10^6,           # oviposition
  delta_E = 1/2.5,    # hatching rate of eggs into larvae
  delta_L = 1/7,   # development rate of larvae into pupae
  delta_P = 1/3,    # development rate of pupae into adult
  tau = 1/2,          # fraction of juviniles becoming adult females
  mu_E = 0.36,      # egg mortality rate
  mu_L = 0.34,      # larva mortality rate
  mu_P = 0.17,      # pupal mortality rate
  mu_A = 0.05,      # adult mortality
  chi = 0.01
)

tf <- 20                     # Final time
simName <- "Mosquito Population Model"   # Model name

# Define initial state vector
x0 <- c(E=10, L=1, P=1, A=1)

# Define state-change matrix
nu <- matrix(
  c(1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, -1),
  nrow = 4,
  byrow = TRUE)

# Define state-change matrix for the model with competition
nu2 <- matrix(
  c(1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, -1, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0),
  nrow = 4,
  byrow = TRUE)

# Define propensity functions
a <- c("alpha*phi*A*(1-(A/K))", "delta_E*E", "delta_L*L", "tau*delta_P*P",
       "mu_E*E", "mu_L*L", "mu_P*P", "mu_A*A")

# Define propensity with competition
a2 <- c("alpha*phi*A*(1-(A/K))", "delta_E*E", "delta_L*L", "tau*delta_P*P",
       "mu_E*E", "mu_L*L", "mu_P*P", "mu_A*A", "chi*L")
# # Run simulations with the Direct method
# out <- ssa(
#   x0 = x0,
#   a = a,
#   nu = nu,
#   parms = parms,
#   tf = tf,
#   method = ssa.d(),
#   simName = simName,
#   verbose = FALSE,
#   consoleInterval = 1
# )
# ssa.plot(out, show.title = TRUE, show.legend = TRUE)

# Plot the times series of the Egg state 
#plot(out$data[,1], out$data[,2], col="red", cex=0.2, pch=19)

###---------------------------------------------------------------------

# Run simulations with the Explict tau-leap method
out_tau.leap <- ssa(
  x0 = x0,
  a = a2,
  nu = nu2,
  parms = parms2,
  tf = tf,
  method = ssa.etl(),
  simName = simName,
  verbose = FALSE,
  consoleInterval = 1) 

ssa.plot(out_tau.leap, show.title = TRUE, show.legend = TRUE)

plot(out_tau.leap$data[,1], out_tau.leap$data[,2], col="red", cex=0.2, pch=19, 
     type = "l", main = "Mosquito Population", xlab = "Time (days)", ylab = "Population")
lines(out_tau.leap$data[,1], out_tau.leap$data[,3], col="green", cex=0.2, pch=19)
lines(out_tau.leap$data[,1], out_tau.leap$data[,4], col="blue", cex=0.2, pch=19)
lines(out_tau.leap$data[,1], out_tau.leap$data[,5], col="violet", cex=0.2, pch=19)
legend(x = "topleft",                                # Position
       legend = c("Egg", "Larva", "Pupa", "Adult"),  # Legend texts
       col = c("red", "green", "blue", "violet"),   # Line colors
       lwd = 2)

# Save the data into a data-frame
data = data.frame(out_tau.leap$data)

# Plot with the ggplot function
ggplot(data, aes(x = t, y = E))+
  geom_line()


# To-d0
# 1. Change the constant parameters to functions and see the outcome of the model
# 2. Make one or two temp-dependent and see the outcome.
# 3. Study to be sure on how to include density-dependence into the model

