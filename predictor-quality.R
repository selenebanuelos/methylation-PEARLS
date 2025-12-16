### Author: Selene Banuelos
### Date: 12/11/2025
### Description: Check for differences in predictions between imputation methods
### and calculate percent of missing predictive CpGs per predictor (want <=10%)

# setup
library(dplyr)
library(tidyr)
options(scipen = 999)

# import data
# generated methylation predictors
predictors <- read.csv('data-processed/methylation-predictors.csv')

# missing CpG per predictor
missing_cpg <- read.csv('data-processed/summary_methscore_CpG.csv')

# calculate differences between predictions generated using each imputation method
pred_diff <- predictors %>%
  select(-c(tissue, contains('resid'))) %>%
  pivot_wider(values_from = -c(SampleID, imp_method),
              names_from = imp_method,
              names_glue = "{.value}{imp_method}"
  ) %>%
  pivot_longer(cols = -SampleID,
               names_to = c('predictor', '.value'),
               names_pattern = ('(.*)(mean|knn)')
  ) %>%
  mutate(difference = mean - knn)
# no differences in any predictors

# calculate the % of predictive CpGs for each predictor ########################
perc_missing <- missing_cpg %>%
  mutate(perc_missing_cpg = 100*(nCpG_missing / nCpG_required))
# 130 missing >10% and 28 missing <=10%