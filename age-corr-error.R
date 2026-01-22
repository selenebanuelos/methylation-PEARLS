### Author: Selene Banuelos
### Date: 1/22/2026
### Description: Visualize correlations between DNAm and chronological age (with
### scatter plots) for predictions generated using PedBE and Horvath2 clocks

# setup
library(dplyr)
library(ggplot2)

# import data
DNAm_age <- read.csv('data-processed/dnam-age-sample-info.csv')

# data wrangling

# data visualization