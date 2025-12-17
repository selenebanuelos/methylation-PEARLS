### Author: Selene Banuelos
### Date: 12/15/2025
### Description: Investigate correlation between chronological age & predicted
### epigenetic ages

# setup
library(dplyr)
library(stringr)
library(corrr)
library(ggplot2)

# import data ##################################################################
# methylation predictor summary
missing_cpg <- read.csv('data-processed/summary_methscore_CpG.csv')

# sample information 
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# methylation predictions
predictions <- read.csv('data-processed/dnam-predictors.csv')

# identify predictors that have <=10% of predictive CpGs missing ###############
quality_predictors <- missing_cpg %>%
  # format predictor name for joining downstream
  mutate(predictor = janitor::make_clean_names(predictor)) %>%
  mutate(perc_missing_cpg = 100*(nCpG_missing / nCpG_required)) %>%
  filter(perc_missing_cpg <= 10)

# data wrangling 
sample_long <- sample %>%
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
               names_pattern = '(.)_(.*)',
               names_to = c('timepoint', '.value')
               ) %>%
  select(c(subjectid, age_baseline, specimenid)) %>%
  mutate(specimenid = as.character(specimenid))

# calculate correlation between predicted age and chronological age ############
# get names of useable predictors (<=10% predictive CpGs missing)
predictors <- quality_predictors$predictor

# data wrangling
ages <- predictions %>%
  mutate(sample_id = str_remove(sample_id,'^0+' # leading zeros
                               ),
         sample_id = str_remove(sample_id,'..$' # last two characters
                               )
  ) %>%
  dplyr::rename(specimenid = sample_id) %>%
  select(-contains(c('resid', 'PC'))) %>%
  full_join(., sample_long, by = 'specimenid') %>%
  filter(imp_method == 'knn') # change to mean if needed after speaking with amy
   
# calculate correlations for blood and buccal samples combined
blood_buccal <- ages %>%
  select(c(age_baseline, contains(predictors))) %>%
  correlate() %>%
  focus(age_baseline) %>%
  dplyr::rename(predictor = term) %>%
  # add predictor type back in
  left_join(., 
            select(quality_predictors, c(predictor, Type)),
            by = 'predictor')

# calculate correlations for blood samples
blood_corr <- ages %>%
  filter(tissue == 'blood') %>%
  select(c(age_baseline, contains(predictors))) %>%
  correlate() %>%
  focus(age_baseline) %>%
  dplyr::rename(predictor = term) %>%
  # add predictor type back in
  left_join(., 
            select(quality_predictors, c(predictor, Type)),
            by = 'predictor')
  
# calculate correlations for buccal samples
buccal_corr <- ages %>%
  filter(tissue == 'buccal') %>%
  select(c(age_baseline, contains(predictors))) %>%
  correlate() %>%
  focus(age_baseline) %>%
  dplyr::rename(predictor = term) %>%
  # add predictor type back in
  left_join(., 
            select(quality_predictors, c(predictor, Type)),
            by = 'predictor')

# plot correlations ############################################################
# age predictors
blood_buccal_plot <- blood_buccal %>%
  filter(str_detect(Type, 'Age')) %>% # keep only age predictors
  mutate(age_baseline = signif(age_baseline, digits = 2)) %>%
  mutate(term = factor(predictor,
                       levels = predictor[order(age_baseline)] # order by corr strength
                       )
         ) %>%
  ggplot(aes(x = predictor,
             y = age_baseline)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = age_baseline),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Correlation with chronological age') +
  xlab('epigenetic age') +
  ggtitle('Buccal and blood samples')

blood_plot <- blood_corr %>%
  filter(str_detect(Type, 'Age')) %>% # keep only age predictors
  mutate(age_baseline = signif(age_baseline, digits = 2)) %>%
  mutate(term = factor(predictor,
                       levels = predictor[order(age_baseline)] # order by corr strength
  )
  ) %>%
  ggplot(aes(x = predictor,
             y = age_baseline)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = age_baseline),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Correlation with chronological age') +
  xlab('epigenetic age') +
  ggtitle('Blood samples')

buccal_plot <- buccal_corr %>%
  filter(str_detect(Type, 'Age')) %>% # keep only age predictors
  mutate(age_baseline = signif(age_baseline, digits = 2)) %>%
  mutate(term = factor(predictor,
                       levels = predictor[order(age_baseline)] # order by corr strength
  )
  ) %>%
  ggplot(aes(x = predictor,
             y = age_baseline)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = age_baseline),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Correlation with chronological age') +
  xlab('epigenetic age') +
  ggtitle('Buccal samples')