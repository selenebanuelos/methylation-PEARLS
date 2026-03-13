### Author: Selene Banuelos
### Date: 3/12/2025
### Description: Compare DNAm age deviation from chronological age across two
### timepoints between participants who experienced no PEARLS and those who
### experienced high PEARLS. Linear mixed effects models will be used to analyze
### this relationship in buccal and blood samples, separately. 

# setup
library(dplyr)
library(lme4)

# import data ##################################################################
# DNAm age data and other sample info
data <- read.csv('data-processed/dnam-age-sample-info.csv')

# demographics data (baseline (T2) household income: income_FPL_100)
demo <- read.csv('data-raw/pearls_dataset_2022-07-08.csv')

# data wrangling ###############################################################
# combine all data needed
clean <- demo %>%
  filter(visitnum == 2) %>% # used only baseline data
  select(pearls_id, income_FPL_100) %>% # get household income data
  right_join(., data, by = 'pearls_id') %>% # combine all data together
  # make PEARLS categories no/high
  mutate(pearls = case_when(aces_baseline == 0 ~ 'no',
                            aces_baseline >= 5 ~ 'high'
                            )
         )

# factor all categorical variables
# timepoint
clean$timepoint <- factor(clean$timepoint,
                          levels = c(2, 5)) # reference: 2
# PEARLS status
clean$pearls <- factor(clean$pearls,
                       levels = c('no', 'high')) # reference: no
# sex
clean$sex <- factor(clean$sex, 
                    levels = c(0, 1),
                    labels = c('female', 'male')) # reference: female

# household income below 100% federal poverty level for household of 4 (<25k)
clean$income_FPL_100 <- factor(clean$income_FPL_100,
                               levels = c(0,1),
                               labels = c('no', 'yes')) # reference: no

# linear mixed effects models ##################################################
# blood samples
blood_ped_be <- lmer(ped_be ~ timepoint*pearls + age_baseline + sex + income_FPL_100 + (1 | pearls_id),
                     data = filter(clean, tissue == 'blood'))

blood_horvath2 <- lmer(horvath2 ~ timepoint*pearls + age_baseline + sex + income_FPL_100 + (1 | pearls_id),
                     data = filter(clean, tissue == 'blood'))

# buccal samples
buccal_ped_be <- lmer(ped_be ~ timepoint*pearls + age_baseline + sex + income_FPL_100 + (1 | pearls_id),
                     data = filter(clean, tissue == 'buccal'))

buccal_horvath2 <- lmer(horvath2 ~ timepoint*pearls + age_baseline + sex + income_FPL_100 + (1 | pearls_id),
                       data = filter(clean, tissue == 'buccal'))