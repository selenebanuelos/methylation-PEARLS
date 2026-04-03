### Author: Selene Banuelos
### Date: 1/22/2026
### Description: Create scatterplots to visualize correlations between DNAm age
### and chronological age for predictions generated using PedBE & Horvath2 clocks
### Calculate median absolute error and add to plots

# setup
library(dplyr)
library(ggplot2)
library(rmcorr) # calculate correlation for repeated measures

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
  mutate(pearls = case_when(aces_baseline == 0 ~ 'no',
                          aces_baseline >= 5 ~ 'high'
                          )
         )

# data visualization
################################################################################
# create scatter plots with regression line to display correlation between
# epigenetic age and chronological age across both timepoints for each tissue
plot <- ages_long %>%
  # only plot data from blood samples
  ggplot(aes(x = age,
             y = dnam_age
             )
         ) +
  
  # scatterplot showing relationship between DNAm age and chrono age
  geom_point(aes(color = pearls)) + # color code points by ACEs status
  
  # add regression line to show correlation
  stat_smooth(method = 'lm',
              formula = y ~ x,
              geom = 'smooth'
              ) +
  
  # stratify plot by epigenetic clock
  facet_grid(vars(clock), vars(tissue)) + 
  
  # add in reference line showing 100% correlation
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'darkgrey') +
  
  # formatting
  labs(title = 'Correlation and median absolute error between epigenetic and chronological age ',
       x = 'Chronological Age',
       y = 'Epigenetic Age') +
  theme_minimal()
plot