### Author: Selene Banuelos
### Date: 5/28/2026
### Description: Summarise DNAm age estimates 

# setup
library(dplyr)

# import data
################################################################################
# DNAm age predictions and sample information
dnam_age <- read.csv('data-processed/dnam-age-sample-info.csv')

# data wrangling
################################################################################
# calculate mean DNAm age for each tissue/timepoint/clock combination
means <- dnam_age %>%
  select(tissue, timepoint, horvath_age, horvath2, ped_be) %>%
  group_by(tissue, timepoint) %>%
  summarise(across(everything(), list(mean = mean)))

# output 
################################################################################
# save means in csv
write.csv(means, 'data-processed/dnam-estimates-summary.csv', row.names = FALSE)