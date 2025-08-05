#Clear workspace
rm(list = ls())

#set working directory
#setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/param-funct-fit")

# Import packages into rcd

library(tidyverse)
library(dplyr)
library(chron)
library(ggplot2)

# Functional forms and parametrization [data taken from Ewing et al, 2016 
#and Beck-Johnson et al, 2013]

# Data on development rates 
egg.data <- data.frame("temp" = c(10, 12, 13, 15, 17, 20, 21, 26, 27, 30, 31),
                       "development_E.sq" = c(0.1, NA, NA, 0.3, 0.35, 0.4, NA, 0.5, NA, 0.8, 1.0),
                       "development_E.ang" = c(NA, 0.2, 0.22, 0.3, 0.32, 0.5, NA, 0.52, 1.0, NA, NA)
                      # "merge" = 0.1, 0.2, 0.22, ...)
)

plot(egg.data$temp, egg.data$development_E.sq, pch=0)
points(egg.data$temp, egg.data$development_E.ang, pch = 17, col = "blue", ylim = c(0, 1))


new.df.1 <- egg.data[,1:2] 
new.df.1$sym <- rep("paper 1", 11)
new.df.2 <- egg.data[,c(1,3)] 
new.df.2$sym <- rep("paper 2", 11)

new.df.1 <- new.df.1 %>% 
  rename(
    development_E.ang = development_E.sq,
  )

new.df.total <- rbind(new.df.1, new.df.2)

#visualize the new data.frame
plot(new.df.total$temp, new.df.total$development_E.ang)
#rbindlist(list(new.df.1,new.df.2))

ggplot(new.df.total, aes(x=temp, y=development_E.ang)) +
  geom_jitter(aes(colour=factor(sym)))+
  xlim(0, 35)

new.df.total<- filter(new.df.total, !is.na(temp), !is.na(development_E.ang)) 

# Fitting a power function to the data using the nls() function
x <- new.df.total$temp
y <- new.df.total$development_E.ang
m <- nls(y ~ a * I(x^b), data = new.df.total, start = list(a = 0.1, b = 0.1))
fit <- fitted(m)

#visualizing with the base plot
#generate range of numbers starting from 0 and ending at 40
xx <- seq(0, 40)
plot(x, y)
lines(x, fit, col = "red", lwd = 2)

#visualizing with ggplot2
ggplot(new.df.total, aes(x=temp, y=development_E.ang)) +
  geom_jitter(aes(colour=factor(sym)))+
  geom_line(aes(x=temp, y=fit, color = "green"))

summary(m)
###########################################################################

# Fitting a polynomial to the egg development dataset
x <- new.df.total$temp
y <- new.df.total$development_E.ang

#first order 
fit.linear <- fit <- lm(y~poly(x, raw = TRUE))
#second degree
fit <- lm(y~poly(x, 2, raw = TRUE))
#third degree
fit1 <- lm(y~poly(x,3,raw=TRUE))
#fourth degree
fit2 <- lm(y~poly(x,4,raw=TRUE))
#generate range of numbers starting from 0 and ending at 40
xx <- seq(0, 40)
plot(x,y,pch=19)
lines(x, egg, col="red")
lines(xx, predict(fit.linear, data.frame(x=xx)), col="red")
lines(xx, predict(fit, data.frame(x=xx)), col="red")
lines(xx, predict(fit1, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="blue")

summary(fit.linear)
summary(fit)
summary(fit1)
summary(fit2)

#visualizing with ggplot2
ggplot(new.df.total, aes(x=temp, y=development_E.ang)) +
  geom_jitter(aes(colour=factor(sym)))+
  geom_line(aes(x=temp, y=predict(fit1), colour = "fit"))
  xlim(0, 35)

 ####################################################################
  new.df.total <- new.df.total[order(new.df.total$development_E.ang),]
  x <- new.df.total$temp
  y <- new.df.total$development_E.ang
  
 
  
  plot(x, y)
  xx <- seq(0, 40, by=0.1)
  a <- 100
  b <- 27
  c <- 1
  y1 <-  c*exp(-(xx-b)^2/a)
  
  plot(xx, y1, type = "l")
  points(x, y)
  
  
  m <- nls(y ~ c*exp(-(x-b)^2/a), data = new.df.total, start = list(a = 100, b = 27, c=1))
  
  plot(x, y)
  lines(x,predict(m),col="red")
  
  summary(m)
  coef(m)
  
  ggplot(new.df.total, aes(x=x, y=y)) +
    geom_jitter()+
    geom_line(aes(x=x, y=predict(m)))
  
  # Just t0 double check the function with the parameter values and data
  a <- unname(coef(m)[1])
  b <- unname(coef(m)[2])
  c <- unname(coef(m)[3])
  
  y2 <-  c*exp(-(x-b)^2/a)
  summary(m)
  
  plot(x, y)
  lines(x, y2, col = "red")
  
  
# Data 
new.df.total
  
# Plot the data
plot(new.df.total$temp, new.df.total$development_E.ang, pch = 17)
lines(x, y2, col = "red")

# This is the function that produces the egg development rate. Use it in running the model. 
   egg.dev <- function(t){
    egg = ifelse(temp(t) <= 10, 0.05, 4.049403 *exp(-(temp(t)-75.098187)^2/1337.666814))
    return(egg)
   }
   
