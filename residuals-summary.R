### Author: Selene Banuelos
### Date: 1/14/2025
### Description: Summarize & plot residuals (age acceleration residuals)

# setup
library(dplyr)
library(ggplot2)

# import data
################################################################################
# generated methylation predictions
predictions <- read.csv('data-processed/methylation-predictions.csv')

# sample information (contains ACEs score)
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# data wrangling
################################################################################
# standardize format of clock names
names(predictions) <- janitor::make_clean_names(names(predictions))

# names of chronological age predictors
chrono_age <- c('horvath_age', 
                'horvath2', 
                'hannum_age', 
                'epigenetic_age_zhang', 
                'c_age', 
                'ped_be'
                )

# clean up sample info for joining with methylation predictions 
sample_long <- sample %>%
  tidyr::pivot_longer(cols = c(T2_specimenid, T5_specimenid),
               names_pattern = '(.*)_(.*)',
               names_to = c('timepoint', '.value')
  ) %>%
  mutate(specimenid = as.character(specimenid)) %>%
  dplyr::rename(pearls_id = subjectid)

# join methylation predictions (& residuals) with sample information
residuals <- predictions %>%
  # create specimenid var to join with sample info downstream
  mutate(specimenid = stringr::str_remove(sample_id,'^0+' # leading zeros
                                          ),
         specimenid = stringr::str_remove(specimenid,'..$' # last two characters
                                          )
         ) %>%
  # keep only predictions for chronological age & sample info
  select(specimenid, age, imp_method, tissue, contains(chrono_age)) %>%
  # keep only residuals of chronological age predictions & sample info
  select(specimenid, age, imp_method, tissue, contains('resid')) %>%
  # add in participant id info to match with demo data
  full_join(., sample_long, by = 'specimenid') %>%
  # keep only predictions generated using KNN for missing CpG imputation
  filter(imp_method == 'knn') %>%
  # create categorical ACEs var (no ACEs/high ACEs) from ACE score at BASELINE
  mutate(aces = case_when(aces_baseline == 0 ~ "no",
                          aces_baseline != 0 ~ "high"
                          )
         )

# reformat residuals data in long form for plotting below
resid_long <- residuals %>%
  tidyr::pivot_longer(
    cols = contains('resid'),
    names_to = c('clock'),
    names_pattern = '(.*)resid',
    values_to = 'resid'
  )

# data visualization
################################################################################
# create 2 plots, one for each tissue type
# create box and whisker plots to show deviation from 0, stratify on ACEs status
# blood plot
resid_long %>%
  filter(clock == 'horvath2' | clock == 'ped_be_') %>%
  filter(tissue == 'blood') %>%
  ggplot(aes(x = aces,# stratify on ACEs status
             y = resid,
             color = clock
             )
         ) +
  geom_boxplot() +
  facet_grid(clock ~ timepoint) +
  ggtitle('Age acceleration residuals in blood samples')
  
# buccal plot
resid_long %>%
  filter(clock == 'horvath2' | clock == 'ped_be_') %>%
  filter(tissue == 'buccal') %>%
  ggplot(aes(x = aces,# stratify on ACEs status
             y = resid,
             color = clock
  )
  ) +
  geom_boxplot() +
  facet_grid(clock ~ timepoint) +
  ggtitle('Age acceleration residuals in buccal samples')