### Author: Selene Banuelos
### Date: 3/12/2025
### Description: Compare DNAm age deviation from chronological age across two
### timepoints between participants who experienced no PEARLS and those who
### experienced high PEARLS. Linear mixed effects models will be used to analyze
### this relationship in buccal and blood samples, separately. 

# setup
library(dplyr)

# import data ##################################################################
# DNAm age data and other sample info
data <- read.csv('data-processed/dnam-age-sample-info.csv')

# demographics data (baseline (T2) household income: income_FPL_100)
demo <- read.csv('data-raw/pearls_dataset_2022-07-08.csv')

# data wrangling ###############################################################
clean <- demo %>%
  filter(visitnum == 2) %>% # used only baseline data
  select(pearls_id, income_FPL_100) %>% # get household income data
  right_join(., data, by = 'pearls_id') %>% # combine all data together
  # make PEARLS categories no/high
  mutate(pearls = case_when(aces_baseline == 0 ~ 'no',
                            aces_baseline >= 5 ~ 'high'
                            )
         )
  
  