# Title: GILLESPIE SSA WITH TIME-DEPENDENT PARAMETERS
# Date: 01/12/2021
# Author: J. Baafi

# Clear workspace
rm(list = ls())

# Set working directory
setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/stochastic-simulations-R")

# Load packages
library(GillespieSSA)
library(ssar)
library(ggplot2)

# GILLESPIE SSA WITH TIME-DEPENDENT PARAMETERS url(https://githubmemory.com/repo/RodrigoZepeda/ssar)
# ------------------------------------------
# Eample 1: First, we run the model with constant parameters
# ------------------------------------------

# Get initial parameters
params <- c(
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


# Initial data must be imputed as a matrix.
X <- matrix(c(E=10, L=1, P=1, A=1), nrow = 1)

#The propensity vector should also be in matrix form:
v <- matrix(
  c(1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, -1),
  nrow = 4,
  byrow = TRUE)

# The propensity scores must also be a matrix-valued function depdendent on 3 parameters: 
# time (t), the state of the system (X) and additional parameters (params) which we 
# discuss later.
pfun <- function(t,X,params){ cbind(params[1]*params[2]*X[, 4]*(1-(X[, 4]/params[3])),
                                    params[4]*X[, 1],
                                    params[5]*X[, 2],
                                    params[7]*params[6]*X[, 3],
                                    params[8]*X[, 1], 
                                    params[9]*X[, 2],
                                    params[10]*X[, 3],
                                    params[11]*X[, 4]) }

# The time for the simulation and number of simulations can be specified:
tmin <- 0
tmax <- 20
nsim <- 2

#Simulate and plot
simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, kthsave = 10, 
                  plot.sim = TRUE,
                  title = "Mosquito Population Model",
                  xlab = "Time", ylab = "Number of individuals")

# With the base plot
plot(simulation$Time, simulation$Var1, "l", col = "red")
lines(simulation$Time, simulation$Var2, col = "blue")
lines(simulation$Time, simulation$Var3, col = "green")
lines(simulation$Time, simulation$Var4, col = "violet")

# -------------------------------------------------------
# Eample 2: Mosquito Population model with time dependent parameters
# ------------------------------------------------------

# ----------------------------------------------------------------------
#EXAMPLE 6
#-----------------------
#Time dependent EA model

# The stage-structured mosquito population model is defined as:
# dE/dt = phi(t)*A*(1-A\K) - delta_E(t)*E
# dA/dt = delta_E(t)*E - delta_A(t)*A
set.seed(123)

#Initial parameters
k <- 100000
params <- c(k = k)
X <- matrix(c(E = 10, A = 1), ncol = 2)
pfun <- function(t, X, params){
  
  #Value to return
  matreturn <- matrix(NA, nrow = length(t), ncol = 3)
  
  #Create birth function
  phi <- function(t){ return(6*exp(t/10)) }

  #Create death function
  deltaE <- function(t){ return(0.3*exp(t/80))}
  
  #Create infectives function
  deltaA <- function(t){ return( 0.3*exp(t/80))}
  
  #Estimate values
  matreturn[,1] <- phi(t)*X[, 2]*(1-X[, 2]/params[1])
  matreturn[,2] <- deltaE(t)*X[, 1]
  matreturn[,4] <- deltaA(t)*X[, 2]
  
  #Return
  return(matreturn)
}

v <- matrix(c(1,-1, 0, 0, 1, -1), nrow = 2, byrow = TRUE)
tmin <- 0
tmax <- 5
nsim <- 2

#Simulate the values
simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim = nsim, print.time = FALSE, 
                  plot.sim = TRUE)

#Plot using ggplot2 library
ggplot(data = simulation, aes(x = Time, y = Var2, group=as.factor(Simulation))) +
  geom_line(aes(color = as.factor(Simulation))) + theme_bw() + 
  theme(legend.position="none") +
  ggtitle(paste0("SIS example; Infected cases ", nsim, " simulations")) + 
  xlab("Time") + ylab("Individuals") +
  geom_vline(xintercept = tmax, linetype = 2)

# try this
set.seed(123)
simulation1 <- ssa(X, pfun, v, params, tmin, tmax, nsim = 10, print.time = FALSE, 
                   plot.sim = FALSE, maxiter = 5000, keep.file = TRUE,
                   fname = "sim1.txt")

set.seed(123)
simulation2 <- ssa(X, pfun, v, params, tmin, tmax, nsim = 10, print.time = FALSE, 
                   plot.sim = FALSE, maxiter = 5000, kthsave = 10, keep.file = TRUE,
                   fname = "sim2.txt")


