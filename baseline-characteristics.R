### Author: Selene Banuelos
### Date: 3/16/2026
### Description: Summarize baseline participant characteristics by PEARLS

# setup
library(dplyr)
library(table1)

# import data
################################################################################
# demographics data
demo <- read.csv('data-raw/pearls_dataset_2022-07-08.csv')

# DNAm age predictions and sample information
dnam_age <- read.csv('data-processed/dnam-age-sample-info.csv')

# data wrangling
################################################################################
# data on household income at T2 (typically 1 month from T2)
income <- demo %>%
  filter(visitnum == 2) %>%
  select(pearls_id,
         income_FPL_100
  )

# combine all data together
characteristics <- dnam_age %>%
  # add in caregiver education at baseline
  left_join(income, by = 'pearls_id') %>%
  # format sex and PEARLS variables
  mutate(pearls = case_when(aces_baseline == 0 ~ 'None',
                            aces_baseline >= 5 ~ 'High'),
         sex = factor(sex,
                      levels = c(0,1),
                      labels = c('Female', 'Male')
         ),
         income_FPL_100 = factor(income_FPL_100,
                                 levels = c(0,1),
                                 labels = c('No', 'Yes')
         )
  ) %>%
  # keep only variables we want in table
  select(pearls_id, sex, age_baseline, income_FPL_100, pearls) %>%
  distinct()

# create labels for variable names to display in table 
label(characteristics$age_baseline) <- 'Age (years)'
label(characteristics$sex) <- 'Sex'
label(characteristics$income_FPL_100) <- 'Household income below 100% FPL (<25k)'
label(characteristics$pearls) <- 'PEARLS'

# data visualization ###########################################################
# table stratified by PEARLS (no/high) status
table1(~ age_baseline + sex + income_FPL_100 | pearls,
       data = characteristics)