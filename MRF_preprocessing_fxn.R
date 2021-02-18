###############################
# Project: Parameter estimation for 3D spatiotempporal Markov Random Field
# Purpose:
# Date: January 21 2021
# Author: Lesley Chapman
# References: MRF.R - Michael Baron
###############################
library(tidyverse)
setwd("./Research_Assistantship/Code")

'
Extract following from dataset:
------------------------------
- Date
- Location : Latitude and Longitude
- omit entries with NA

input
-----
csv file listing incident reports: date, time, coordinates

output
------
dataframe listing: date, latitude, longitude
'

extract_data <- function(file) {
  df <- readr::read_csv(file)
  # Show number of entries
  # nrow(df)
  
  #Extract: Date, Lat, Long
  df %>% 
    select("callDateTime", "location") %>%
    separate(callDateTime, c('Date','Time'), "\\s+") %>%
    separate("location", c('Latitude', 'Longitude'), sep=",") %>%
    mutate(Latitude = gsub("\\(", "", Latitude)) %>% 
    mutate(Longitude = gsub("\\)", "", Longitude)) %>%
    select (-c(Time)) %>%
    mutate(Latitude = as.numeric(Latitude)) %>%
    mutate(Longitude = as.numeric(Longitude)) %>%
    drop_na() -> df
  
  # Recode Date, Longitude, Latitude
  df <- transform(df,Date2 = as.numeric(factor(Date)))
  # df <- transform(df,Latitude2 = as.numeric(factor(Latitude)))
  # df <- transform(df,Longitude2 = as.numeric(factor(Longitude)))

  return(df)
}

# return_coor <- function(df, lat_1, lat_2, lon_1, lon_2) {
#   df %>% 
#     mutate(Latitude = as.numeric(Latitude)) %>%
#     mutate_at("Latitude", funs(round(., 2))) %>%
#     mutate(Longitude = as.numeric(Longitude)) %>%
#     mutate_at("Longitude", funs(round(., 2))) %>%
#     select(Latitude,Longitude,Latitude2,Longitude2)-> df
#   
#   df <- df[!duplicated(df), ]
#   
#   head(df)
# 
#   # lat <- ifelse(grepl(lat_1 == df$Latitude, df$Latitude2, lat_1))
#   #               
#   # lat
# }

'
Generate delta values

input
-----
df : date, lat, long
h: number specifying IQR coefficient
x_cntr, y_cntr: specify coordinates for city center; ex: xx = c(-76.67,-76.57); yy = c(39.27,39.32)

output
------
Delta X
Delta Y
'

find_deltaValues <- function(file, h, x_cntr, y_cntr){
  # Count number of events per day
  file %>% 
    group_by(Date) %>%
    mutate(call_ct = n()) %>% 
    mutate(callgap_ct = ifelse(call_ct < 10, 
                               "Gap in timeline!!! Only",
                               "entries on day")) -> df
  
  # ID outliers outside of Q1 and Q3 IQR
  y1 = quantile(df$Longitude,.25) - h*IQR(df$Longitude)
  y2 = quantile(df$Longitude,.75) + h*IQR(df$Longitude)
  x1 = quantile(df$Latitude,.25) - h*IQR(df$Latitude)
  x2 = quantile(df$Latitude,.75) + h*IQR(df$Latitude)

  # Generate Delta Values
  Nbins = 30
  DeltaX = (x_cntr[2]-x_cntr[1])/(Nbins +1)
  DeltaY = (y_cntr[2]-y_cntr[1])/(Nbins +1)

  delta_values <- list(DX = DeltaX, DY = DeltaY)
  return(delta_values)
}



#######################################################################
### Estimate MRF parameters - 2D king size neighborhood ###############
#######################################################################
 
find_counts <- function(df, Nbins, x_cntr, y_cntr, DeltaX, DeltaY){
  # Xbin and Ybin
  # Capture the square in which each event has occured
  # The following represents squares in the lattice
  # Need to account for every square in which an event has occured
  Xbin <- ceiling((df$Latitude - x_cntr[1])/DeltaX)
  Ybin <- ceiling((df$Longitude - y_cntr[1])/DeltaY)
  Tbin <- df$Date2 - min(df$Date2) + 1
  Ndays <- max(Tbin)
  
  g <- cbind(Xbin, Ybin, Tbin)
  g <- as.data.frame(g)
  return(head(g))

  Nlocations = Nbins*Nbins
  Events = rep(0, Ndays*Nlocations)
  dim(Events) = c( Ndays, Nbins, Nbins )
  IndexEvents = which( Xbin >= 1 & Xbin <= Nbins & Ybin >= 1 & Ybin <= Nbins )
  
  #Transform matrix into X and Y vectors
  LocationX = rep(0,Nlocations)
  LocationY = LocationX

  for (i in 1:Nbins){
    LocationX[(Nbins*(i-1)+1) : (Nbins*i)] = i
    LocationY[(Nbins*(i-1)+1) : (Nbins*i)] = seq(1,Nbins)
  }

  # Neighboring locations
  i0 = seq(2,(Nbins-1));   # The point itself
  i1 = seq(1,(Nbins-2));   # Its left or bottom neighbor
  i2 = seq(3,Nbins);       # Its right or top neighbor


  # Define the 3D array:
  for (i in IndexEvents){
    Events[ Tbin[i], Xbin[i], Ybin[i] ] = 1
  }

  # Consider each day, save the logistic regression coefficients
  Parameters = matrix(rep(0, Ndays*10),Ndays,10)

  for (t in 1:Ndays){
    DayEvents = Events[t,,];       # Nbins*Nbins matrix of 0s and 1s
    # Define 10 Delta-statistics for each location,
    # Each of them is a (Nbins-2)*(Nbins-2) matrix of 0s and 1s
    
    # Using each location as is...
    N0  = DayEvents[i0,i0]
    N1 = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i0,i1]*(df$Date2+1))
    N2     = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i0,i2]*(df$Date2+1))
    N3     = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i1,i0]*(df$Date2+1))
    N4    = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i2,i0]*(df$Date2+1))
    N5    = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i1,i0] + DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i1,i1])
    N6   = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i2,i0] + DayEvents[i0,i2]*DayEvents[i2,i2] + DayEvents[i1,i0]*DayEvents[i1,i1])
    N7  = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i1,i0] + DayEvents[i0,i2]*DayEvents[i2,i1] + DayEvents[i0,i2]*DayEvents[i1,i2])
    N8   = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0] + DayEvents[i0,i1]*DayEvents[i2,i1] + DayEvents[i1,i0]*DayEvents[i1,i2])
    N9  = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i2,i0]*DayEvents[i2,i1] +
                                  DayEvents[i0,i1]*DayEvents[i1,i0]*DayEvents[i1,i1] + DayEvents[i0,i2]*DayEvents[i1,i0]*DayEvents[i1,i2])
    
    # Replacing each location with the opposite...
    n0  = (1-DayEvents[i0,i0])
    n1 = (1-DayEvents[i0,i0])*(DayEvents[i2,i0] + DayEvents[i1,i0])
    n2     = (1-DayEvents[i0,i0])*(DayEvents[i0,i2] + DayEvents[i0,i1])
    n3     = (1-DayEvents[i0,i0])*(DayEvents[i2,i2] + DayEvents[i1,i1])
    n4    = (1-DayEvents[i0,i0])*(DayEvents[i2,i1] + DayEvents[i1,i2])
    n5    = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i1,i0] + DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i1,i1])
    n6   = (1-DayEvents[i0,i0])*(DayEvents[i0,i1]*DayEvents[i2,i0] + DayEvents[i0,i2]*DayEvents[i2,i2] + DayEvents[i1,i0]*DayEvents[i1,i1])
    n7  = (1-DayEvents[i0,i0])*(DayEvents[i0,i1]*DayEvents[i1,i0] + DayEvents[i0,i2]*DayEvents[i2,i1] + DayEvents[i0,i2]*DayEvents[i1,i2])
    n8   = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i2,i0] + DayEvents[i0,i1]*DayEvents[i2,i1] + DayEvents[i1,i0]*DayEvents[i1,i2])
    n9  = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i2,i0]*DayEvents[i2,i1] +
                                      DayEvents[i0,i1]*DayEvents[i1,i0]*DayEvents[i1,i1] + DayEvents[i0,i2]*DayEvents[i1,i0]*DayEvents[i1,i2])

    # And taking differences = Delta-statistics...
    D0  = -as.numeric(as.vector(Nalpha - nalpha))
    D1 = -as.numeric(as.vector(Nlambda - nlambda))
    D2     = -as.numeric(as.vector(Nmu - nmu))
    D3     = -as.numeric(as.vector(Npi - npi))
    D4    = -as.numeric(as.vector(Nrho - nrho))
    D5    = -as.numeric(as.vector(Neta - neta))
    D6   = -as.numeric(as.vector(Nzeta - nzeta))
    D7  = -as.numeric(as.vector(Ntheta - ntheta))
    D8   = -as.numeric(as.vector(Niota - niota))
    D9  = -as.numeric(as.vector(Nkappa - nkappa))

    EventResponse = as.factor(Nalpha);

    logreg = glm( EventResponse ~ D0 + D1 + D2
                  + D3 + D4 + D5 + D6 + D7 + D8 + D9,
                  family="binomial",
                  control=list(maxit = 500) )

    Parameters[t,] = coef(logreg)[2:11]
  }

}



