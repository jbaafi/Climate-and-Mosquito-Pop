# This code solves a stage-structured ODE model for mosquito population dynamics.

# clears the workspace
rm(list = ls())

#Load R package for solving Differential equations
library(deSolve)
library(chron)
library(tidyverse)
library(dplyr)
library(ggplot2)
# Set working directory for the script
setwd("/Users/jbaafi/Documents/Research@MUN/model parameterization data")

climate.df <- read.csv("climate.df.csv", header = TRUE)

#Define a function for temperature
t <- climate.df$Days.Since.Origin
Temp <- climate.df$Mean.Temp
#a <- mean(climate.df$Mean.Temp)
a <- 8.863282
b1 <- -8.7635
b2 <- b1
#t <- seq(1, 365)

Temp.func <- function(t){
Temp <- a + b1*sin(2*pi*t/365) + b2*cos(2*pi*t/365)
return(Temp)
}

plot(t, Temp.func(t))

#Define the parameters
#Recruitment into the egg state. Data taken from Ewing et al, 2020.
eggs <- function(Temp.func){
  egg <- 0.003903*Temp.func(t)^2 -0.156807*Temp.func(t)+1.619391
  return(egg)
}


plot(t, Temp.func(t))
plot(Temp.func(t), eggs(Temp.func(t)))



