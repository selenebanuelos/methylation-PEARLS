### Author: Selene Banuelos
### Date: 12/15/2025
### Description: Investigate correlation between chronological & predicted ages

# setup
library(dplyr)
library(stringr)

# import data
# methylation predictor summary
missing_cpg <- read.csv('data-processed/summary_methscore_CpG.csv')

# sample information 
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# methylation predictions
predictions <- read.csv('data-processed/methylation-predictors.csv')

# select predictors that have <=10% of predictive CpGs missing
quality_predictors <- missing_cpg %>%
  filter(100*(nCpG_missing / nCpG_required) <= 10) %>%
  # format predictor name for joining downstream
  mutate(predictor = str_replace(predictor, ' ', '.'),
         predictor = str_replace(predictor, '-', '.'),
         predictor = str_extract(predictor, 
                                          '^[^ ]+' #keep everything until 1st space
                                          )
         )

# data wrangling 
sample_long <- sample %>%
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
               names_pattern = '(.)_(.*)',
               names_to = c('timepoint', '.value')
               ) %>%
  select(c(subjectid, age_baseline, specimenid)) %>%
  mutate(specimenid = as.character(specimenid))


# plot correlation between predicted age and chronological age
# get names of useable predictors
predictors <- quality_predictors$predictor

ages <- predictions %>%
  mutate(SampleID = str_remove(SampleID,'^0+' # leading zeros
                               ),
         SampleID = str_remove(SampleID,'..$' # last two characters
                               )
  ) %>%
  dplyr::rename(specimenid = SampleID) %>%
  full_join(., sample_long, by = 'specimenid') %>%
  select(-contains('resid')) %>%
  select(contains(predictors)) %>%
  select(-contains('PC'))

# plot correlations
         