###############################
# Project: Parameter estimation for 3D spatiotempporal Markov Random Field
# Purpose:
# Date: January 21 2021
# Author: Lesley Chapman
# References: MRF.R - Michael Baron
###############################
library(tidyverse)
setwd("/Volumes/Lesley_Chapman/American_University/Research_Assistantship/Code")

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
    # mutate(Latitude = as.numeric(Latitude)) %>%
    # mutate(Longitude = as.numeric(Longitude)) %>%
    drop_na() -> df
  
  # Recode Date, Longitude, Latitude
  df <- transform(df,Date2 = as.numeric(factor(Date)))
  df <- transform(df,Latitude2 = as.numeric(factor(Latitude)))
  df <- transform(df,Longitude2 = as.numeric(factor(Longitude)))

  return(df)
}

return_coor <- function(df, lat_1, lat_2, lon_1, lon_2) {
  df %>% 
    mutate(Latitude = as.numeric(Latitude)) %>%
    mutate_at("Latitude", funs(round(., 2))) %>%
    mutate(Longitude = as.numeric(Longitude)) %>%
    mutate_at("Longitude", funs(round(., 2))) %>%
    select(Latitude,Longitude,Latitude2,Longitude2)-> df
  
  df <- df[!duplicated(df), ]
  
  head(df)

  # lat <- ifelse(grepl(lat_1 == df$Latitude, df$Latitude2, lat_1))
  #               
  # lat
}

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
  y1 = quantile(df$Longitude2,.25) - h*IQR(df$Longitude2)
  y2 = quantile(df$Longitude2,.75) + h*IQR(df$Longitude2)
  x1 = quantile(df$Latitude2,.25) - h*IQR(df$Latitude2)
  x2 = quantile(df$Latitude2,.75) + h*IQR(df$Latitude2)

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
  # Find standardized locations Xbin and Ybin
  Xbin <- ceiling((df$Latitude2 - x_cntr[1])/DeltaX)
  Ybin <- ceiling((df$Longitude2 - y_cntr[1])/DeltaY)
  Tbin <- df$Date2 - min(df$Date2) + 1
  Ndays <- max(Tbin)

  
  Nlocations = Nbins*Nbins
  Events = rep(0, Ndays*Nlocations)
  dim(Events) = c( Ndays, Nbins, Nbins )
  IndexEvents = which( Xbin >= 1 & Xbin <= Nbins & Ybin >= 1 & Ybin <= Nbins );

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
    Nalpha  = DayEvents[i0,i0]
    Nlambda = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i0,i1]*(df$Date2+1))
    Nmu     = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i0,i2]*(df$Date2+1))
    Npi     = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i1,i0]*(df$Date2+1))
    Nrho    = (DayEvents[i0,i0]*df$Date2)*(DayEvents[i2,i0]*(df$Date2+1))
    Neta    = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i1,i0] + DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i1,i1])
    Nzeta   = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i2,i0] + DayEvents[i0,i2]*DayEvents[i2,i2] + DayEvents[i1,i0]*DayEvents[i1,i1])
    Ntheta  = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i1,i0] + DayEvents[i0,i2]*DayEvents[i2,i1] + DayEvents[i0,i2]*DayEvents[i1,i2])
    Niota   = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0] + DayEvents[i0,i1]*DayEvents[i2,i1] + DayEvents[i1,i0]*DayEvents[i1,i2])
    Nkappa  = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i2,i0]*DayEvents[i2,i1] +
                                  DayEvents[i0,i1]*DayEvents[i1,i0]*DayEvents[i1,i1] + DayEvents[i0,i2]*DayEvents[i1,i0]*DayEvents[i1,i2])

    # Replacing each location with the opposite...
    nalpha  = (1-DayEvents[i0,i0])
    nlambda = (1-DayEvents[i0,i0])*(DayEvents[i2,i0] + DayEvents[i1,i0])
    nmu     = (1-DayEvents[i0,i0])*(DayEvents[i0,i2] + DayEvents[i0,i1])
    npi     = (1-DayEvents[i0,i0])*(DayEvents[i2,i2] + DayEvents[i1,i1])
    nrho    = (1-DayEvents[i0,i0])*(DayEvents[i2,i1] + DayEvents[i1,i2])
    neta    = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i1,i0] + DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i1,i1])
    nzeta   = (1-DayEvents[i0,i0])*(DayEvents[i0,i1]*DayEvents[i2,i0] + DayEvents[i0,i2]*DayEvents[i2,i2] + DayEvents[i1,i0]*DayEvents[i1,i1])
    ntheta  = (1-DayEvents[i0,i0])*(DayEvents[i0,i1]*DayEvents[i1,i0] + DayEvents[i0,i2]*DayEvents[i2,i1] + DayEvents[i0,i2]*DayEvents[i1,i2])
    niota   = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i2,i0] + DayEvents[i0,i1]*DayEvents[i2,i1] + DayEvents[i1,i0]*DayEvents[i1,i2])
    nkappa  = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i2,i0]*DayEvents[i2,i1] +
                                      DayEvents[i0,i1]*DayEvents[i1,i0]*DayEvents[i1,i1] + DayEvents[i0,i2]*DayEvents[i1,i0]*DayEvents[i1,i2])

    # And taking differences = Delta-statistics...
    Dalpha  = -as.numeric(as.vector(Nalpha - nalpha))
    Dlambda = -as.numeric(as.vector(Nlambda - nlambda))
    Dmu     = -as.numeric(as.vector(Nmu - nmu))
    Dpi     = -as.numeric(as.vector(Npi - npi))
    Drho    = -as.numeric(as.vector(Nrho - nrho))
    Deta    = -as.numeric(as.vector(Neta - neta))
    Dzeta   = -as.numeric(as.vector(Nzeta - nzeta))
    Dtheta  = -as.numeric(as.vector(Ntheta - ntheta))
    Diota   = -as.numeric(as.vector(Niota - niota))
    Dkappa  = -as.numeric(as.vector(Nkappa - nkappa))

    EventResponse = as.factor(Nalpha);

    logreg = glm( EventResponse ~ Dalpha + Dlambda + Dmu
                  + Dpi + Drho + Deta + Dzeta + Dtheta + Diota + Dkappa,
                  family="binomial",
                  control=list(maxit = 500) )

    Parameters[t,] = coef(logreg)[2:11]
  }

}



