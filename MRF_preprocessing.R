###############################
# Project: Parameter estimation for 3D spatiotempporal Markov Random Field
# Purpose:
# Date: January 21 2021
# Author: Lesley Chapman
# References: MRF.R - Michael Baron
###############################
source('/Volumes/Lesley_Chapman/American_University/Research_Assistantship/Code/MRF/MRF_preprocessing_fxn.R')
library(tidyverse)
setwd("/Volumes/Lesley_Chapman/American_University/Research_Assistantship/Code")
 
# Variables
Nbins <- 30

# Files
file <- "./data/Baltimore_911_Calls_for_Service.csv"

# Run code
df_event <- extract_data(file)
delta_values <- find_deltaValues(df_event, Nbins, c(-76.67,-76.57), c(39.27,39.32))
deltaX <- delta_values$DX
deltaY <- delta_values$DY

# Estimate MRF parameters - 2D king size neighborhood
x <- find_counts(df_event, Nbins, c(-76.67,-76.57), c(39.27,39.32), deltaX, deltaY)
x
