# Title:  Stochastic compartmental model for mosquito population dynamics
# Author: Joseph Baafi
# Date:   March 24, 2022.

# Clear workspace
rm(list = ls())

# Set working directory
# setwd("directory here")

# Load packages 
pacman::p_load(pacman, GillespieSSA, ssar, ggplot2)

# The stage-structured mosquito population model is defined as:
# dE/dt = b*rhoA*A - (delta_E + mu_E)*E
# dL/dt = delta_E*E - (delta_L + mu_L)*L
# dP/dt = delta_L*L - (delta_P + mu_P)*P
# dA/dt = tau*(delta_P*P) - mu_A*A

#Get initial parameters
params <- c(b = 100, rhoA = 3, delta_E = 0.5, delta_L = 0.14,
            delta_P = 0.5, mu_E = 0.56, mu_L = 0.44, mu_P = 0.37, mu_A = 0.4)

X <- matrix(c(100, 100, 100, 100), ncol = 4)

#Propensity function
pfun       <- function(t, X, params){ cbind(params[1]*params[2]*X[,4], 
                                            params[3]*X[,1], 
                                            params[4]*X[,2],
                                            params[5]*X[,3],
                                            params[6]*X[,1],
                                            params[7]*X[,2], 
                                            params[8]*X[,3],
                                            params[9]*X[,4]) }
#The propensity vector should also be in matrix form:
v <- matrix(
  c(1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, -1),
  nrow = 4,
  byrow = TRUE)

#Simulate
simulation <- ssa(X, pfun, v, params, nsim = 2, tmax = 2,
                  title = "Model with constant parameters",
                  xlab = "Time", ylab = "Number of individuals")



# I will get the above code running fine and start replacing the parameters with 
# functions. 




set.seed(322)
X          <- matrix(c(N=500), nrow = 1)   #Initial values
v          <- matrix( c(+1, -1), ncol = 2) 

f <- function(t){sin(t*pi)}
pfun       <- function(t,X,params){ cbind(2 * X[, 1], 
                                          (2 + f(t)*X[, 1]/1000)*X[, 1]) }

simulation <- ssa(X, pfun, v, tmin = 2, tmax = 10, nsim = 2, 
                  title = "Time-dependent Logistic Growth: Example2", 
                  xlab = "Time", ylab = "Individuals")


# Lotka-Volterra model
#Set seed
set.seed(3289650)

#Get initial parameters
params <- c(a = 3, b = 0.01, c = 2)
X <- matrix(c(100, 100), ncol = 2)

#Propensity function
pfun       <- function(t, X, params){ cbind(params[1]*t*X[,1] + 1, 
                                            params[2]*X[,1]*X[,2], 
                                            params[3]*X[,2]) }
#Propensity score
v          <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)

#Simulate
simulation <- ssa(X, pfun, v, params, nsim = 2, tmin = 0, tmax = 10,
                  title = "Example 4: Time-dependent Lotka-Volterra",
                  xlab = "Time", ylab = "Number of individuals")

sim        <- read.table("My_simulation.txt",  header = TRUE)

ggplot(data = sim, aes(x = Time, group = as.factor(Simulation))) +
  geom_line(aes(y = Var1, color = "Prey")) +
  geom_line(aes(y = Var2, color = "Predator")) +
  ggtitle("Example 4: Lotka Volterra with ggplot2") + 
  xlab("Time") + ylab("Individuals") +
  scale_color_manual("Creature", 
                     values = c("Prey" = "deepskyblue4","Predator" = "tomato3"))

# SIS model
#Initial parameters
k          <-  24576.5529836797
delta      <-  0.0591113454895868 + 0.208953907151055
gamma_ct   <-  0.391237630231631
params     <- c(k = k, delta = delta, gamma_ct = gamma_ct)
X          <- matrix(c(S = 1000000, I = 1000), ncol = 2)
pfun       <- function(t, X, params){
  
  #Value to return
  matreturn  <- matrix(NA, nrow = length(t), ncol = 6)
  
  #Create birth function
  lambda     <- function(t){ return(4.328e-4 - (2.538e-7)*t - 
                                      (3.189e-7)*sin(2 * t * pi/52) - 
                                      (3.812e-7)*cos(2 * t * pi/52) ) }
  
  #Create death function
  mu         <- function(t){ return(9.683e-5 + (1.828e-8)*t + 
                                      (2.095e-6)*sin(2 * t * pi/52) - 
                                      (8.749e-6)*cos(2 * t * pi/52))}
  
  #Create infectives function
  beta_fun   <- function(t){ return( 0.479120824267286 + 
                                       0.423263042762498*sin(-2.82494252560096 + 2*t*pi/52) )}
  
  #Estimate values
  matreturn[,1] <- lambda(t)*(X[,1] + X[,2])
  matreturn[,2] <- mu(t)*X[,1]
  matreturn[,3] <- beta_fun(t)*X[,1]*X[,2]/(1 + params[1]*X[,2])
  matreturn[,4] <- mu(t)*X[,2]
  matreturn[,5] <- params[2]*X[,2]
  matreturn[,6] <- params[3]*X[,2]
  
  #Return
  return(matreturn)
  
}
v          <- matrix(c(1,-1, -1, 0, 0, 1, 0, 0, 1, -1, -1, -1), nrow = 2, byrow = TRUE)
tmin       <- 0
tmax       <- 10
nsim       <- 100

simulation1 <- ssa(X, pfun, v, params, tmin, tmax, nsim = 2, print.time = FALSE, 
                   plot.sim = TRUE, keep.file = TRUE,
                   fname = "sim1.txt")

plot(simulation1$Time, simulation1$Var1, "l")
plot(simulation1$Time, simulation1$Var2, "l")

# PLots
plot(t, 4.328e-4 - (2.538e-7)*t - (3.189e-7)*sin(2 * t * pi/52) - 
       (3.812e-7)*cos(2 * t * pi/52), "l", main = "Birth function")

plot(t, 9.683e-5 + (1.828e-8)*t + (2.095e-6)*sin(2 * t * pi/52) - 
       (8.749e-6)*cos(2 * t * pi/52), main = "Death function", "l")

plot(t, 0.479120824267286 + 
       0.423263042762498*sin(-2.82494252560096 + 2*t*pi/52), main = "Infectives Rate function", "l")




