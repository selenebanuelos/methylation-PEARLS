### Author: Selene Banuelos
### Date: 3/16/2026
### Description: Summarize baseline participant characteristics by PEARLS status

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
# data on household income at T2 (typically 1 month from T1)
income <- demo %>%
  filter(visitnum == 2) %>%
  select(pearls_id, income_FPL_100)

# data on caregiver education at T1
edu <- demo %>%
  filter(visitnum == 1) %>%
  select(pearls_id, caregiver_edu_4groups)

# combine all data together
characteristics <- dnam_age %>%
  # add in houshold income at baseline
  left_join(income, by = 'pearls_id') %>%
  # caregiver education at baseline
  left_join(edu, by = 'pearls_id') %>%
  # format sex, PEARLS, and household income variables
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
  # keep only variables wanted in table
  select(pearls_id, 
         sex, 
         age_baseline, 
         income_FPL_100, 
         pearls, 
         caregiver_edu_4groups) %>%
  distinct()

# create labels for variable names to display in table 
label(characteristics$age_baseline) <- 'Age (years)'
label(characteristics$sex) <- 'Sex'
label(characteristics$income_FPL_100) <- 'Household income below 100% FPL (<25k)'
label(characteristics$caregiver_edu_4groups) <- 'Caregiver highest educational attainment'
label(characteristics$pearls) <- 'PEARLS'

# data visualization ###########################################################
# check distribution of chronological age
hist(characteristics$age_baseline) # overall
filter(characteristics, pearls == 'None') %>% pull(age_baseline) %>% hist()
filter(characteristics, pearls == 'High') %>% pull(age_baseline) %>% hist()
# looks like there's a right skew (more younger children)

# table stratified by PEARLS (no/high) status
table1(~ age_baseline + sex + income_FPL_100 + caregiver_edu_4groups | pearls,
       data = characteristics,
       render.continuous = c(.='Median [Min, Max]'))