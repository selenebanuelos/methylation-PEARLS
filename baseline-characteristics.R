### Author: Selene Banuelos
### Date: 3/16/2025
### Description: Summarize baseline participant characteristics by PEARLS
### Note: only include samples that passed methylation data QC

# setup
library(dplyr)
library(table1)

# import data
################################################################################
# sample information
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# demographics data
demo <- read.csv('data-raw/pearls_dataset_2022-07-08.csv')

# samples that passed qc
passed_qc <- read.csv('data-processed/samples-passed-qc.csv') %>%
  select(specimenid, passed_qc, tissue) %>%
  mutate(specimenid = as.character(specimenid))

# data wrangling
################################################################################
# data on caregiver education at T1 (baseline)
edu <- demo %>%
  filter(visitnum == 1) %>%
  select(pearls_id,
         caregiver_edu_binary, 
         caregiver_edu_4groups
  )

# data on household income at T2 (typically 1 month from T2)
income <- demo %>%
  filter(visitnum == 2) %>%
  select(pearls_id,
         visitnum,
         income_FPL_100
  )

# combine all data together
characteristics <- sample %>%
  # format sex and PEARLS variables
  mutate(pearls = case_when(aces_baseline == 0 ~ 'None',
                            aces_baseline >= 5 ~ 'High'),
         sex = factor(sex,
                         levels = c(0,1),
                         labels = c('Female', 'Male')
                         ),
  ) %>%
  # rename participant id for downstream joining
  rename(pearls_id = subjectid) %>%
  # add in car
         