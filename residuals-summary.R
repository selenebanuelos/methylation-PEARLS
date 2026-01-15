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

# data wrangling
################################################################################
# join methylation predictions (& residuals) with sample information
residuals <- predictions %>%
  mutate(specimenid = str_remove(sample_id,'^0+' # leading zeros
  ),
  specimenid = str_remove(specimenid,'..$' # last two characters
  )
  ) %>%
  select(specimenid, age, imp_method, tissue, contains(chrono_age)) %>%
  select(specimenid, age, imp_method, tissue, contains('resid')) %>%
  # add in participant id info to match with demo data
  full_join(., sample_long, by = 'specimenid') %>%
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
  facet_wrap(~clock, nrow = 1)
  
# buccal plot