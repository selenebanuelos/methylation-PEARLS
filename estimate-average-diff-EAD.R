### Author: Selene Banuelos
### Date: 1/10/2025
### Description: Linear regression analysis to estimate difference in average 
### EAD between those who reported experiencing no PEARLS items to those who 
### reported ≥ 5 PEARLS will be presented within each timepoint.

# setup
library(dplyr)

# import data
################################################################################
# residuals from regressing epigenetic age ~ chronological age for each timepoint
residuals <- read.csv('data-processed/epi-chrono-age-residuals.csv')

# demographics data
demo <- read.csv('data-raw/pearls_dataset_2022-07-08.csv') 

# data wrangling
################################################################################
all_data <- demo %>%
  filter(visitnum == 2) %>%
  select(pearls_id, income_FPL_100) %>%
  right_join(.,
             residuals,
             by = "pearls_id")
buccal_2 <- all_data

# regression
################################################################################
# EAD from PedBE clock
buccal_2 

buccal_5 

blood_2

blood_5

# EAD from Skin & Blood clock
buccal_2

buccal_5

blood_2

blood_5
