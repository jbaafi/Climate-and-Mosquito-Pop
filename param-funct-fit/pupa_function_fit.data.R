#Clear workspace
rm(list = ls())

#set working directory
setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/param-funct-fit")

# Import packages into rcd

library(tidyverse)
library(dplyr)
library(chron)
library(ggplot2)

# Functional forms and parametrization [data taken from Ewing et al, 2016 
#and Beck-Johnson et al, 2013]

# Data on development rates primarily from Ewing et al.
pupa.data <- data.frame("temp" = c(14, 15, 16, 19, 20, 23, 26, 27, 28, 30, 35),
                        "develop_L.circ" = c(NA, 0.17, NA, NA, 0.4, NA, 0.42, 0.6, NA, 0.55, 0.54),
                        "develop_L.ang" = c(0.06, NA, 0.05, 0.2, 0.1, 0.12, NA, NA, 0.25, NA, NA),
                        "develp_L.star" = c(NA, 0.05, NA, NA, 0.16, NA, 0.24, NA, 0.25, NA, NA))

pupa.df.1 <- pupa.data[,1:2] 
pupa.df.1$sym <- rep("paper 1", 11)
pupa.df.2 <- pupa.data[,c(1,3)] 
pupa.df.2$sym <- rep("paper 2", 11)
pupa.df.3 <- pupa.data[,c(1,4)]
pupa.df.3$sym <- rep("paper 3", 11)

pupa.df.1 <- pupa.df.1 %>% 
  rename(
    development.rate = develop_L.circ,
  )

pupa.df.2 <- pupa.df.2 %>% 
  rename(
    development.rate = develop_L.ang,
  )

pupa.df.3 <- pupa.df.3 %>% 
  rename(
    development.rate = develp_L.star,
  )

pupa.df.total <- rbind(pupa.df.1, pupa.df.2, pupa.df.3)

#visualize the new data.frame
plot(pupa.df.total$temp, pupa.df.total$development.rate)
#rbindlist(list(new.df.1,new.df.2))

ggplot(pupa.df.total, aes(x=temp, y=development.rate)) +
  geom_jitter(aes(colour=factor(sym)))

pupa.df.total<- filter(pupa.df.total, !is.na(temp), !is.na(development.rate)) 

# Fitting a polynomial to the larva development dataset
x <- pupa.df.total$temp
y <-pupa.df.total$development.rate
#third degree
fit1 <- lm(y~poly(x,3,raw=TRUE))
#fourth degree
fit2 <- lm(y~poly(x,4,raw=TRUE))
#generate range of numbers starting from 0 and ending at 40
xx <- seq(0, 40)
plot(x,y,pch=19)
lines(xx, predict(fit1, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="blue")

summary(fit1)
summary(fit2)

#visualizing with ggplot2
ggplot(pupa.df.total, aes(x=temp, y=development.rate)) +
  geom_jitter(aes(colour=factor(sym)))+
  geom_line(aes(x=temp, y=predict(fit1), colour = "fit"))+
  geom_line(aes(x=temp, y=predict(fit2), colour = "fit2"))


#Fitting a bell-shaped function to the data using the nls() function (non-linear least squares)
pupa.df.total <- pupa.df.total[order(pupa.df.total$temp),]
#pupa.df.total <- pupa.df.total[order(pupa.df.total$development.rate),]

x <- pupa.df.total$temp
y <-pupa.df.total$development.rate
xx <- seq(0, 40, by=0.1)
a <- 60
b <- 27
c <- 0.6
y1 <-  c*exp(-(x-b)^2/a)

plot(x, y1, type = "l") 
points(x, y)

m <- nls(y ~ c*exp(-(x-b)^2/a), data = pupa.df.total, start = list(a = 60, b = 27, c=0.6))

plot(x, y)
lines(x, predict(m))

coef(m)

pupa <- function(x){
  p <- 0.5920232 *exp(-(x-39.9070020)^2/339.6702830)
  return(p)
  }

plot(x, y)             
lines(x, pupa(x))


# The function to be used to run the model is as follows: I used pupa2 just to confirm that it is working fine
pupa2 <- function(x){
  p <- ifelse(x <= 10 | x >= 40, 0.03, 0.5920232 *exp(-(x-39.9070020)^2/339.6702830))
  return(p)
}

lines(x, pupa2(x), col = "red")
