#Clear workspace
rm(list = ls())

#set working directory
setwd("/Users/jbaafi/Documents/Research@MUN/model parameterization data")


library(tidyverse)
library(dplyr)
library(chron)
library(ggplot2)

# Functional forms and parametrization [data taken from Ewing et al, 2016 
#and Beck-Johnson et al, 2013]

# Data on adult death rates  primarily taken from Beck-Johnson et al,
temp <- c(5, 10, 15, 20, 25, 30, 35, 40)
mortality.paper1 <- c(0.358, 0.190, 0.053, 0.042, 0.042, 0.114, 0.232, 1.283)
mortality.paper2 <- c(0.358, 0.144, 0.053, 0.049, 0.053, 0.107, 0.215, 0.849)

adult.mortality.df <- data.frame(temp, mortality.paper1, mortality.paper2)

#visualize data with base plot function
plot(adult.mortality.df$temp, adult.mortality.df$mortality.paper1, pch = 5)
points(adult.mortality.df$temp, adult.mortality.df$mortality.paper2, pch = 8, col = "red")

#visualize data with ggplot2
ggplot(data=adult.mortality.df)+
  geom_point(mapping = aes(x=temp, y=mortality.paper1, color = "red"), shape = 5)+
  geom_point(mapping = aes(x=temp, y=mortality.paper2, color = "blue"), shape=8)

# Fitting a model to the data
adult.df.1 <- adult.mortality.df[,1:2]
adult.df.1$sym <- rep("paper 1", 8)
adult.df.2 <- adult.mortality.df[,c(1,3)] 
adult.df.2$sym <- rep("paper 2", 8)

adult.df.1 <- adult.df.1 %>% 
  rename(
    mortality.rate = mortality.paper1,
  )

adult.df.2 <- adult.df.2 %>% 
  rename(
    mortality.rate = mortality.paper2,
  )

adult.df.total <- rbind(adult.df.1, adult.df.2)

#visualize the new adult mortality data.frame
plot(adult.df.total$temp, adult.df.total$mortality.rate)
#rbindlist(list(new.df.1,new.df.2))

ggplot(data = adult.df.total, mapping = aes(x=temp, y=mortality.rate))+
  geom_jitter(aes(colour=factor(sym)))


# Fitting a polynomial to the adult mortality dataset
x <- adult.df.total$temp
y <- adult.df.total$mortality.rate
#third degree
fit1 <- lm(y~poly(x,3,raw=TRUE))
#fourth degree
fit2 <- lm(y~poly(x,4,raw=TRUE))
#generate range of numbers starting from 0 and ending at 40
xx <- seq(0, 50, by=0.1)
plot(x,y,pch=19)
lines(xx, predict(fit1, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="blue")

summary(fit1)
summary(fit2)

#visualizing with ggplot2
ggplot(adult.df.total, aes(x=temp, y=mortality.rate)) +
  geom_jitter(aes(colour=factor(sym)))+
  geom_line(aes(x=temp, y=predict(fit1), colour = "fit"))+
  geom_line(aes(x=temp, y=predict(fit2), colour = "fit2"))

################################################################################
#Fitting a non-linear function to the data using the nls() function
adult.df.total <- adult.df.total[order(adult.df.total$temp),]

x <- adult.df.total$temp
y <- adult.df.total$mortality.rate

plot(x, y)

a <- 0.045
b <- 21
c <- 0.0035
y1 <-  c*(x-b)^2+a

c_0 <- 0.0886
c_1 <- 21.211
c_2 <- 14.852
  
y2 <- c_0*exp((x-c_1)/c_2)^4

plot(x, y2, "l", col = "red")
lines(x, y1) 
points(x, y)


m <- nls(y ~ c*(x-b)^2+a, data = adult.df.total, start = list(a = 0.045, b = 21, c=0.0035))

m2 <- nls(y ~ c_A*exp(((x-T_A)/d_A)^4), data = adult.df.total, 
          start = list(c_A = 0.0886, T_A = 21, d_A = 14.852))


plot(x, y)
lines(x, predict(m))
lines(x, predict(m2), col = "red")

func.data <- data.frame(x, predict(m2))

summary(func.data)

summary(m)
summary(m2)
###############################################################################
#Model parameters
a.temp <- 8.9231 
b1.temp <- -8.7635
b2.temp <- b1.temp

t <-  seq(0, 365)

# Defining temperature function
temp <- function(t){
  temp = a.temp + b1.temp*sin(2 * pi * t/365) + b2.temp*cos(2 * pi * t/365)
  return(temp)
}

data.temp <- data.frame(t, temp(t))

#Plot of temperature as a function of time (365 days)
ggplot(data = data.temp, mapping = aes(x=t, y=temp(t)))+
  geom_line()

c_A <- 0.08841
T_A <- 21.24746
d_A <- 14.92552


#look at this
f <- function(t){
  func = ifelse(temp(t) <= 5, 0.36003054, c_A*exp(((temp(t)-T_A)/d_A)^4))
  return(func)
}

data.f <- data.frame(t, temp(t), f(t))

plot(data.f$t, data.f$f.t., type = "l")

summary(data.f)


