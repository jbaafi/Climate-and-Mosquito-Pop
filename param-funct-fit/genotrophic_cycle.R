#Clear workspace
rm(list = ls())

# Load data
pacman::p_load(tidyverse, dplyr, chron, ggplot2)

# Data on genotrophic cycle rates  primarily taken from Beck-Johnson et al, and Lardeux et al, 2008
temp <- c(15, 20, 23, 25, 27, 29, 31, 33, 35)
genotrophic.cycle <- c(0.0781, 0.1563, 0.1961, 0.2857, 0.3030, 0.3704, 0.4762, 0.4167, 0.4348)

# Put it in a dataframe
df <- data.frame(temp=temp, gen.cycle=genotrophic.cycle)

#visualize data with base plot function
plot(df$temp, df$gen.cycle, pch=16, col="darkblue", cex=1)

temp <- df$temp

gen.func <- function(alpha, beta1, beta2){
  alpha*exp(-((temp-beta1)^2)/beta2)
}

model1 <- nls(df$gen.cycle ~ gen.func(alpha, beta1, beta2), data = df, start = list(alpha=0.46, beta1 = 31, beta2 = 100))

summary(model1)
coef(model1)

phi <-  0.443111*exp(-((temp - 34.388928)^2)/185.611107) # Exponential function for genotrophic cycle

plot(df$temp, df$gen.cycle, pch=16, col="darkblue", cex=1, ylim = c(0, 0.5))
lines(df$temp, phi, col = "red")


# Fit power function to data same data
power <- function(alpha, beta) {alpha*temp^beta }

model2 <- nls(df$gen.cycle ~ power(alpha, beta), data = df, start = list(alpha = 0.001, beta = 1.7))

summary(model2)
coef(model2)

phi2 <- 0.0007419766*temp^1.8246405203 # Power function for genotrophic cycle
# Plot
plot(df$temp, df$gen.cycle, pch=16, col="darkblue", cex=1)
lines(df$temp, phi2, "l")

AIC(model1)
AIC(model2)
BIC(model1)
BIC(model2)


# Estimating the parameters with MLE of the function alpha*exp(-(temp-beta1)^2/beta2)
temp <- df$temp

LL <- function(alpha, beta1, beta2, mu, sigma) {
  R = df$gen.cycle - alpha*exp(-((temp-beta1)^2)/beta2)
  #
  R = suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
  #
  -sum(R)
}

library(bbmle)

fit <- mle2(LL, start = list(alpha = 0.3, beta1 = 34, beta2 = 100, mu = 0, sigma = 1))

summary(fit)
AIC(fit)

phi2 <- coef(fit)[1]*exp(-((temp-coef(fit)[2])^2)/coef(fit)[3])
plot(df$temp, df$gen.cycle, pch=16, col="darkblue", cex=1, ylim = c(0, 0.5))
lines(df$temp, phi2, "l")
# This seem not to fit very well
