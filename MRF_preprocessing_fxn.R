###############################
# Project: Parameter estimation for 3D spatiotempporal Markov Random Field
# Purpose:
# Date: January 21 2021
# Author: Lesley Chapman
# References: MRF.R - Michael Baron
###############################
library(tidyverse)
setwd("/Volumes/Lesley_Chapman/American_University/Research_Assistantship/Code/3D_MarkovRandomField")

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
  
  ## Recode Date, Longitude, Latitude
  df <- transform(df,Date_adj = as.numeric(factor(Date)))
  # df <- transform(df,Latitude2 = as.numeric(factor(Latitude)))
  # df <- transform(df,Longitude2 = as.numeric(factor(Longitude)))
  
  return(df)
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
### Estimate MRF parameters - 3D neighborhood ###############
#######################################################################
'
Find all possible combinations of Xbin, Ybin and Tbin
'

generate_bins <- function(df, Nbins){
  # Import dataframe and recode: date, latitude, and longitude
  df %>%
    #Sort Dates in ascending order
    mutate(Date = as.Date(Date, "%m/%d/%Y")) %>%
    arrange(Date) -> city_df

  city_df$Latitude <- ceiling(city_df$Latitude)
  city_df$Longitude <- ceiling(city_df$Longitude)
  
  #Recode date, latitude, and longitude
  city_df <- transform(city_df,Date_adj = as.numeric(factor(Date)))
  city_df <- transform(city_df,Lat = as.numeric(factor(Latitude)))
  city_df <- transform(city_df,Lon = as.numeric(factor(Longitude)))
  city_df %>% select(-c(Date, Latitude, Longitude)) -> city_df
  
  # Create dataframe with all combinations of Xbin, Ybin, Tbin
  Ndays <- max(city_df$Date_adj)
  N = Ndays*Nbins*Nbins
  N
  
  T = rep(0,N)
  X=T
  Y=T
  data.frame(T,X,Y)
  
  for (t in 1:Ndays){
    T[((t-1)*Nbins*Nbins + 1) : (t*Nbins*Nbins)]  = t
    for (x in 1:Nbins) {
      X[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = x #Give all Values for X
      Y[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = 1:Nbins #Give all Values for Y
    }
  }
  
  null_df <- data.frame(T,X,Y)
  
  null_df %>%
    rename(Date_adj = T) %>%
    left_join(city_df, by=c('Date_adj')) %>%
    mutate(Z = ifelse(X == Lat & Y == Lon, 1, 0)) %>%
    rename(Tbin = Date_adj) %>%
    rename(Xbin = X) %>%
    rename(Ybin = Y) %>%
    select(-c(Lat, Lon))-> x
  
  return(head(x))
  
}

xytz_bin_test <- function(n, bin){
  Ndays = n # Fill this in based on computation in MRF code
  Nbins = bin
  
  N = Ndays*Nbins*Nbins
  N
  
  T = rep(0,N)
  X=T
  Y=T
  data.frame(T,X,Y)
  
  for (t in 1:Ndays){
    T[((t-1)*Nbins*Nbins + 1) : (t*Nbins*Nbins)]  = t
    for (x in 1:Nbins) {
      X[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = x #Give all Values for X
      Y[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = y #Give all Values for Y
    }
  }
  
  Z = rep(0,N)
  Z[(t-1)*Nbins^2 + (x-1)*Nbins + y] = 1
  return(head(data.frame(T,X,Y,Z)))
}

find_counts <- function(df, Nbins, x_cntr, y_cntr, DeltaX, DeltaY){
  # Xbin and Ybin
  # Capture the square in which each event has occured
  # The following represents squares in the lattice
  # Need to account for every square in which an event has occured
  Xbin <- ceiling((df$Longitude - x_cntr[1])/DeltaX)
  Ybin <- ceiling((df$Latitude - y_cntr[1])/DeltaY)
  Tbin <- df$Date_adj - min(df$Date_adj) + 1
  Ndays <- max(Tbin)
  
  Xbin_adj <- data.frame(transform(as.numeric(factor(Xbin))))
  Ybin_adj <- data.frame(transform(as.numeric(factor(Ybin))))
  df_bin <- cbind.data.frame(Xbin_adj,Ybin_adj, Tbin)
  df_bin <- df_bin %>% 
    add_column(Z=0) 
  df_bin <- df_bin %>% rename(Xbin = X_data)
  df_bin <- df_bin %>% rename(Ybin = X_data.1)
   return(head(df_bin))
  
  
  
  Nlocations = Nbins*Nbins
  Z = rep(0, Ndays*Nlocations)
  dim(Z) = c( Ndays, Nbins, Nbins )
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
    Z[ Tbin[i], Xbin[i], Ybin[i] ] = 1
  }
  
  # Consider each day, save the logistic regression coefficients
  Parameters = matrix(rep(0, Ndays*10),Ndays,10)
  
  for (t in 1:Ndays){
    DayEvents = Z[t,,];       # Nbins*Nbins matrix of 0s and 1s
    # Define 10 Delta-statistics for each location,
    # Each of them is a (Nbins-2)*(Nbins-2) matrix of 0s and 1s
    
    # Using each location as is...
    N0  = Z[i0,i0,i0] # Take sum here?
    N1 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(1 + Nbins*Nbins + Nbins) : N]) #North East
    N2 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(1 + Nbins*Nbins + 1) : N]) #North
    N3 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(1 + Nbins*Nbins + 1 + Nbins) : N]) #North West
    N4 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(1 + Nbins*Nbins - Nbins) : N]) #West
    N5 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(1 + Nbins*Nbins) : N]) #Same location
    N7 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(Nbins*Nbins + Nbins) : (N-1)]) #East
    N7 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(Nbins*Nbins - Nbins) : (N-1)]) #Southwest
    N8 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(Nbins*Nbins) : (N-1)]) #South
    N9 = sum(Z[1:(N-Nbins*Nbins-Nbins)] * Z[(Nbins*Nbins + Nbins) : (N-1)]) #Southeast 
    # N1 = (DayEvents[i0,i0]*df$Date_adj)*(DayEvents[i0,i1]*(df$Date_adj+1))
    # N2     = (DayEvents[i0,i0]*df$Date_adj)*(DayEvents[i0,i2]*(df$Date_adj+1))
    # N3     = (DayEvents[i0,i0]*df$Date_adj)*(DayEvents[i1,i0]*(df$Date_adj+1))
    # N4    = (DayEvents[i0,i0]*df$Date_adj)*(DayEvents[i2,i0]*(df$Date_adj+1))
    #N5    = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i1,i0] + DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i1,i1])
    # N6   = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i2,i0] + DayEvents[i0,i2]*DayEvents[i2,i2] + DayEvents[i1,i0]*DayEvents[i1,i1])
    # N7  = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i1,i0] + DayEvents[i0,i2]*DayEvents[i2,i1] + DayEvents[i0,i2]*DayEvents[i1,i2])
    # N8   = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0] + DayEvents[i0,i1]*DayEvents[i2,i1] + DayEvents[i1,i0]*DayEvents[i1,i2])
    #N9  = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i2,i0]*DayEvents[i2,i1] +
    #            DayEvents[i0,i1]*DayEvents[i1,i0]*DayEvents[i1,i1] + DayEvents[i0,i2]*DayEvents[i1,i0]*DayEvents[i1,i2])
    
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
