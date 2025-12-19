### Author: Selene Banuelos
### Date: 12/15/2025
### Description: Select which clocks to use based on criteria
### Investigate correlation between chronological age & epigenetic ages

# setup
library(dplyr)
library(stringr)
library(corrr)
library(ggplot2)

# import data ##################################################################
# information on required CpGs per predictor
pred_info <- read.csv('data-processed/summary_methscore_CpG.csv')

# generated methylation predictions
predictions <- read.csv('data-processed/dnam-predictors.csv')

# participant information (including chronological age)
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# calculate % of missing predictive CpGs and filter for clocks only (no protein)
clocks <- pred_info %>%
  # format predictor name for joining downstream
  mutate(predictor = janitor::make_clean_names(predictor),
         # calculate % of missing CpGs
         perc_missing_cpg = 100*(nCpG_missing / nCpG_required)
  ) %>%
  # filter for clocks/epigenetic age predictors
  filter(Type == 'DNAm Age') %>%
  # remove principal component (PC) versions of clocks since they are redundant
  filter(!grepl('based on PC component', Description))
# warning is from weird character in some descriptions

# get vector of clock names
clock_names <- clocks$predictor

# clean up sample data
sample_long <- sample %>%
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
               names_pattern = '(.)_(.*)',
               names_to = c('timepoint', '.value')
               ) %>%
  select(c(subjectid, age_baseline, specimenid)) %>%
  mutate(specimenid = as.character(specimenid))

# clean up prediction data & merge to participant information
ages <- predictions %>%
  mutate(sample_id = str_remove(sample_id,'^0+' # leading zeros
                               ),
         sample_id = str_remove(sample_id,'..$' # last two characters
                               )
  ) %>%
  dplyr::rename(specimenid = sample_id) %>%
  select(-contains(c('resid', 'PC'))) %>%
  full_join(., sample_long, by = 'specimenid') %>%
  filter(imp_method == 'knn')

# calculate correlations for blood samples
blood_corr <- ages %>%
  filter(tissue == 'blood') %>%
  select(c(age_baseline, contains(clock_names))) %>%
  correlate() %>%
  focus(age_baseline) %>%
  dplyr::rename(predictor = term,
                corr_chrono_age = age_baseline) %>%
  # add in % of missing predictive CpGs
  mutate(perc_missing_cpg = clocks$perc_missing_cpg,
         cat_perc_missing_cpg = case_when(
           perc_missing_cpg <= 10 ~ '% <= 10',
           perc_missing_cpg > 10 & perc_missing_cpg <= 20 ~ '10 < % <= 20',
           perc_missing_cpg > 20 & perc_missing_cpg <= 30 ~ '20 < % <= 30',
           perc_missing_cpg > 30 ~ '% > 30'
         )
  )
  
# calculate correlations for buccal samples
buccal_corr <- ages %>%
  filter(tissue == 'buccal') %>%
  select(c(age_baseline, contains(clock_names))) %>%
  correlate() %>%
  focus(age_baseline) %>%
  dplyr::rename(predictor = term,
                corr_chrono_age = age_baseline) %>%
  # add in % of missing predictive CpGs
  mutate(perc_missing_cpg = clocks$perc_missing_cpg,
         cat_perc_missing_cpg = case_when(
           perc_missing_cpg <= 10 ~ '% <= 10',
           perc_missing_cpg > 10 & perc_missing_cpg <= 20 ~ '10 < % <= 20',
           perc_missing_cpg > 20 & perc_missing_cpg <= 30 ~ '20 < % <= 30',
           perc_missing_cpg > 30 ~ '% > 30'
           )
         )

# data visualization ###########################################################
# order bar plots by % missing cpg (low to high) and put above bar?
# bar plots showing correlation between epigenetic & chrono age in blood
blood_plot <- blood_corr %>%
  # format numbers
  mutate(corr_chrono_age = signif(corr_chrono_age, digits = 2),
         perc_missing_cpg = signif(perc_missing_cpg, digits = 2)
         ) %>%
  ggplot(aes(
    # reorder bars in ascending order based on % missing predictive CpGs
    x = reorder(predictor, perc_missing_cpg), 
    y = corr_chrono_age
    )) +
  geom_bar(stat = 'identity') +
  # label correlations between ages within bars
  geom_text(aes(label = corr_chrono_age),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Correlation with chronological age') +
  xlab('epigenetic age') +
  ggtitle('Blood samples')

# bar plots showing correlation between epigenetic & chrono age in buccal swabs
buccal_plot <- buccal_corr %>%
  mutate(corr_chrono_age = signif(corr_chrono_age, digits = 2)) %>%
  ggplot(aes(
    # bars in ascending order based on % missing predictive CpGs
    x = reorder(predictor, perc_missing_cpg), 
    y = corr_chrono_age)
    ) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = corr_chrono_age),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Correlation with chronological age') +
  xlab('epigenetic age') +
  ggtitle('Buccal swab samples')

# create bar plots showing mean absolute error & order by % missing CpGs