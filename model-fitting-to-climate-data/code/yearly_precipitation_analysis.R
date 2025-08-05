# Title: Analyze of climate data and fit functions to temperature and precipitation data
# Author: Joseph Baafi

#Clear workspace
rm(list = ls())

#set working directory
setwd("/Users/jbaafi/Documents/climate-and-mosquitoes/model-fitting-to-climate-data/code")

# Import packages into r
#library(tidyverse)
library(dplyr)
library(chron)
library(ggplot2)
library(gridExtra)
library(grid)


# Load 2011-2016 climate data
climate_11 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2011_P1D.csv", header = TRUE)
climate_12 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2012_P1D.csv", header = TRUE)
climate_13 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2013_P1D.csv", header = TRUE)
climate_14 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2014_P1D.csv", header = TRUE)
climate_15 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2015_P1D.csv", header = TRUE)
climate_16 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2016_P1D.csv", header = TRUE)

# Use select() function from dplyr
select <- dplyr::select

# Select needed columns
df_11 <- climate_11 %>% 
  select(Station.Name, Year:Day, Total.Precip..mm.)

# Formate time data in a form that chron can understand
daily_dates <- dates(paste(df_11$Month, 
                           df_11$Day, 
                           df_11$Year, sep="/")) 
df_11$chron.date <- chron(dates=daily_dates, 
                              origin. = c(month = 1,day = 1,
                                          year = 2011)) # setting date into chron
df_11 <- arrange(df_11, chron.date) # ordering dates and adding it to data-frame

start.date <- as.Date("2011-01-01") # The starting date
df_11$days.since.origin <- (as.numeric(as.Date(df_11$chron.date) - start.date)) # Producing a sequence from smallest to largest

# Count NANs column-wise
sapply(df_11, function(x) sum(is.na(x)))

# Format data to exclude the NANs 
df_11<- filter(df_11, !is.na(Total.Precip..mm.)) 

# Plot data
plot1 <- ggplot(data = df_11, mapping = aes(x=chron.date, y=Total.Precip..mm.))+
  geom_line(colour = "steelblue")+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Precipitation")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# --------------------------------------------------

# Subset 2012 data
df_12 <- climate_12 %>% 
  select(Station.Name, Year:Day, Total.Precip..mm.)

# Format time data in a form that chron can understand 
daily_dates.12 <- dates(paste(df_12$Month, 
                           df_12$Day, 
                           df_12$Year, sep="/")) 
df_12$chron.date <- chron(dates=daily_dates.12, 
                          origin. = c(month = 1,day = 1,
                                      year = 2012)) # set date into chron
df_12 <- arrange(df_12, chron.date) 

start.date.12 <- as.Date("2012-01-01") # The starting date
df_12$days.since.origin <- (as.numeric(as.Date(df_12$chron.date) - start.date.12)) # Producing a sequence from smallest to largest

# Formate data to exclude the NANs 
df_12<- filter(df_12, !is.na(Total.Precip..mm.)) 

# Time series plot of temperature 
plot2 <- ggplot(data = df_12, mapping = aes(x=chron.date, y=Total.Precip..mm.))+
  geom_line(colour = "steelblue")+
  ylim(0, 60)+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Precipitation")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# --------------------------------------------------

# 2013 data
df_13 <- climate_13 %>% 
  select(Station.Name, Year:Day, Total.Precip..mm.)

# Format time data
daily_dates.13 <- dates(paste(df_13$Month, 
                              df_13$Day, 
                              df_13$Year, sep="/")) 
df_13$chron.date <- chron(dates=daily_dates.13, 
                          origin. = c(month = 1,day = 1,
                                      year = 2013)) # setting date into chron
df_13 <- arrange(df_13, chron.date) # ordering dates and adding it to data-frame

start.date.13 <- as.Date("2013-01-01") # The starting date
df_13$days.since.origin <- (as.numeric(as.Date(df_13$chron.date) - start.date.13)) # Producing a sequence from smallest to largest

# Delete NANs
df_13<- filter(df_13, !is.na(Total.Precip..mm.)) 

# Plot of 2013 data 
plot3 <- ggplot(data = df_13, mapping = aes(x=chron.date, y=Total.Precip..mm.))+
  geom_line(colour = "steelblue")+
  ylim(0, 60)+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Precipitation")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# --------------------------------------------------

# Plots in a grid form
#grid.arrange(plot1, plot2, plot3, ncol=1)

# 2014 data
df_14 <- climate_14 %>% 
  select(Station.Name, Year:Day, Total.Precip..mm.)

# Formate time data in a form that chron can understand 
daily_dates.14 <- dates(paste(df_14$Month, 
                              df_14$Day, 
                              df_14$Year, sep="/")) 
df_14$chron.date <- chron(dates=daily_dates.14, 
                          origin. = c(month = 1,day = 1,
                                      year = 2014)) # setting date into chron
df_14 <- arrange(df_14, chron.date) # ordering dates and adding it to data-frame

start.date.14 <- as.Date("2014-01-01") # The starting date
df_14$days.since.origin <- (as.numeric(as.Date(df_14$chron.date) - start.date.14)) # Producing a sequence from smallest to largest


# Exclude the NANs 
df_14<- filter(df_14, !is.na(Total.Precip..mm.)) 

# Time series plot of 2014 data
plot4 <- ggplot(data = df_14, mapping = aes(x=chron.date, y=Total.Precip..mm.))+
  geom_line(colour = "steelblue")+
  ylim(0, 60)+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Precipitation")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# --------------------------------------------------

grid.arrange(plot1, plot2, plot3, plot4, ncol=1)


# 2015 data
df_15 <- climate_15 %>% 
  select(Station.Name, Year:Day, Total.Precip..mm.)

# Formatting time data in a form that chron can understand 
daily_dates.15 <- dates(paste(df_15$Month, 
                              df_15$Day, 
                              df_15$Year, sep="/")) 
df_15$chron.date <- chron(dates=daily_dates.15, 
                          origin. = c(month = 1,day = 1,
                                      year = 2015)) # setting date into chron
df_15 <- arrange(df_15, chron.date) # ordering dates and adding it to data-frame

start.date.15 <- as.Date("2015-01-01") # The starting date
df_15$days.since.origin <- (as.numeric(as.Date(df_15$chron.date) - start.date.15)) # Producing a sequence from smallest to largest

# Formating data to exclude the NANs 
df_15 <- filter(df_15, !is.na(Total.Precip..mm.)) 

# Time series plot 
plot5 <- ggplot(data = df_15, mapping = aes(x=chron.date, y=Total.Precip..mm.))+
  geom_line(colour = "steelblue")+
  ylim(0, 60)+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Precipitation")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# --------------------------------------------------

# 2016
df_16 <- climate_16 %>% 
  select(Station.Name, Year:Day, Total.Precip..mm.)

# Formatting time data in a form that chron can understand 
daily_dates.16 <- dates(paste(df_16$Month, 
                              df_16$Day, 
                              df_16$Year, sep="/")) 
df_16$chron.date <- chron(dates=daily_dates.16, 
                          origin. = c(month = 1,day = 1,
                                      year = 2016)) # setting date into chron
df_16 <- arrange(df_16, chron.date) # ordering dates and adding it to data-frame

start.date.16 <- as.Date("2016-01-01") # The starting date
df_16$days.since.origin <- (as.numeric(as.Date(df_16$chron.date) - start.date.16)) # Producing a sequence from smallest to largest

# Formating data to exclude the NANs 
df_16 <- filter(df_16, !is.na(Total.Precip..mm.)) 

# Time series plot
plot6 <- ggplot(data = df_16, mapping = aes(x=chron.date, y=Total.Precip..mm.))+
  geom_line(colour = "steelblue")+
  ylim(0, 60)+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Precipitation")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
# --------------------------------------------------

# Grid plots (ggplot)
grid.arrange(plot1, plot2, plot3, ncol=1, top="Daily precipitation for 2011 - 2013")
grid.arrange(plot4, plot5, plot6, ncol = 1, top="Daily precipitation for 2014 - 2016")


# Combine the above data-frames into one data-frame.
df <- rbind(df_11, df_12, df_13, df_14, df_15, df_16) 

#------------------
# It looks like there's no particular trend in the data in terms of seasonality. 
# I should be able to figure out what period seems to be the wet and the dry season.
# Maybe try plotting data from Africa if that can be easily seen.
# ----------------

# To-do 
# Super-impose one data on the other by using different colours to indicate different years
plot(df_11$days.since.origin, df_11$Total.Precip..mm., type = "l", col = "red")
lines(df_12$days.since.origin, df_12$Total.Precip..mm., col = "blue")
lines(df_13$days.since.origin, df_13$Total.Precip..mm., col = "green")
legend("topleft", legend=c("2011", "2012", "2013"),
      col=c("red", "blue", "green"), lty=1, cex=0.8)

# I should probably forget about the trend in seasons and use one model to represent the data
# Incorporating rainfall into the model:
# Rainfall/precipitation should influence the model and the amount of the rainfall will
# be drawn from the fitted distribution, so you find rainfall functions on the parameters 
# as it is done in Abdelrazak but rainfall will be drawn from these pdf stochasticaly.
