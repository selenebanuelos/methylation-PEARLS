### Author: Selene Banuelos
### Date: 1/22/2026
### Description: Visualize correlations between DNAm and chronological age (with
### scatter plots) for predictions generated using PedBE and Horvath2 clocks

# setup
library(dplyr)
library(ggplot2)

# import data
################################################################################
DNAm_age <- read.csv('data-processed/dnam-age-sample-info.csv')

# data wrangling
################################################################################
# clean up DNAm and chrono age data for plotting
ages_long <- DNAm_age %>%
  select(horvath2, ped_be, age, tissue, aces_baseline, timepoint) %>%
  # make data longer for plotting downstream
  tidyr::pivot_longer(
    cols = c(horvath2, ped_be),
    names_to = 'clock',
    values_to = 'dnam_age'
  ) %>%
  # create categorical ACEs variable
  mutate(aces = case_when(aces_baseline == 0 ~ 'no',
                          aces_baseline >= 5 ~ 'high'
                          )
         )

# data visualization
################################################################################
# create scatter plots with regression line to display correlation between
# epigenetic age and chronological age across both timepoints for each tissue

# blood samples
ages_long %>%
  # only plot data from blood samples
  filter(tissue == 'blood') %>%
  ggplot(aes(x = age,
             y = dnam_age
             )
         ) +
  # scatterplot showing relationship between DNAm age and chrono age
  geom_point(
    aes(color = aces) # color code points by ACEs status
    ) +
  # add regression line to show correlation
  stat_smooth(method = 'lm',
              formula = y ~ x,
              geom = 'smooth'
              ) +
  # stratify plot by epigenetic clock
  facet_wrap(~clock)

# buccal samples
ages_long %>%
  # only plot data from blood samples
  filter(tissue == 'buccal') %>%
  ggplot(aes(x = age,
             y = dnam_age
  )
  ) +
  # scatterplot showing relationship between DNAm age and chrono age
  geom_point(
    aes(color = aces) # color code points by ACEs status
  ) +
  # add regression line to show correlation
  stat_smooth(method = 'lm',
              formula = y ~ x,
              geom = 'smooth'
  ) +
  # stratify plot by epigenetic clock
  facet_wrap(~clock)