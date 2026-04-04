### Author: Selene Banuelos
### Date: 1/22/2026
### Description: Create scatterplots to visualize correlations between DNAm age
### and chronological age for predictions generated using PedBE & Horvath2 clocks
### Calculate median absolute error and add to plots
### Correlations and MAE caclulated for each timepoint separately

# setup
library(dplyr)
library(ggplot2)
library(Metrics) # calculatge meadian absolute error

# import data
################################################################################
DNAm_age <- read.csv('data-processed/dnam-age-sample-info.csv')

# data wrangling
################################################################################
# clean up DNAm and chrono age data for plotting
ages_long <- DNAm_age %>%
  select(pearls_id, horvath2, ped_be, age, tissue, aces_baseline, timepoint) %>%
  # make data longer for plotting downstream
  tidyr::pivot_longer(
    cols = c(horvath2, ped_be),
    names_to = 'clock',
    values_to = 'dnam_age'
  ) %>%
  # create categorical PEARLS variable
  mutate(Pearls = case_when(aces_baseline == 0 ~ 'No',
                            aces_baseline >= 5 ~ 'High'
  ),
  # format clock names for plot
  clock = case_when(clock == 'horvath2' ~ 'Skin & Blood',
                    clock == 'ped_be' ~ 'PedBE'),
  # format tissue names for plot
  tissue = case_when(tissue == 'blood' ~ 'Blood',
                     tissue == 'buccal' ~ 'Buccal')
  )

# calculate DNAm age ~ chrono age correlations 
################################################################################
cor <- ages_long %>%
  # calculate correlations at each timepoint separately for each tissue/clock group
  group_by(timepoint, tissue, clock) %>%
  summarize(corr = round(cor(age, dnam_age),
                         digits = 2
                         )
            )

# calculate DNAm age - chrono age median absolute error 
################################################################################
mae <- ages_long %>%
  # calculate MAE within each timepoint separately for each tissue/clock group
  group_by(timepoint, tissue, clock) %>%
  summarize(mae = round(mdae(actual = age, predicted = dnam_age),
                        digits = 2)
            ) %>%
  ungroup()

# data visualization
################################################################################
# create scatter plots with regression line to display correlation between
# epigenetic age and chronological age across both timepoints for each tissue
baseline <- ages_long %>%
  # only plot data from baseline
  filter(timepoint == 2) %>%
  ggplot(aes(x = age, y = dnam_age)) +
  
  # scatterplot showing relationship between DNAm age and chrono age
  geom_point(shape = 16, aes(color = Pearls)) + # color code points by ACEs status
  
  # add regression line to show correlation
  stat_smooth(method = 'lm', formula = y ~ x) +
  
  # stratify plot by epigenetic clock
  facet_grid(clock ~ tissue) +
  
  # add in reference line showing 100% correlation
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'darkgrey') +
  
  # add in correlation coefficients (r)
  geom_text(data = filter(cor, timepoint == 2),
            aes(label = paste0('r = ', corr), 
                x = 2,
                y = 16)) +
  
  # add in median absolute errors
  geom_text(data = filter(mae, timepoint == 2), 
            aes(label = paste0('MAE = ', mae), 
                x = 2.4, 
                y = 14)) +
  
  # formatting
  labs(title = 'Correlation and median absolute error (MAE) between epigenetic and chronological age ',
       x = 'Chronological Age',
       y = 'Epigenetic Age') +
  
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        plot.caption = element_text(hjust = 0),
        legend.position = 'bottom') +
  
  labs(caption = paste('Figure. Performance of epigenetic clocks from samples at baseline.',
                       '\n',
                       'Dashed line is reference for perfect linear relationship',
                       '\n',
                       'r = Pearson correlation coefficient')
       )
baseline

followup <- ages_long %>%
  # only plot data from baseline
  filter(timepoint == 5) %>%
  ggplot(aes(x = age, y = dnam_age)) +
  
  # scatterplot showing relationship between DNAm age and chrono age
  geom_point(shape = 16, aes(color = Pearls)) + # color code points by ACEs status
  
  # add regression line to show correlation
  stat_smooth(method = 'lm', formula = y ~ x) +
  
  # stratify plot by epigenetic clock
  facet_grid(clock ~ tissue) +
  
  # add in reference line showing 100% correlation
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'darkgrey') +
  
  # add in correlation coefficients (r)
  geom_text(data = filter(cor, timepoint == 5),
            aes(label = paste0('r = ', corr), 
                x = 6.5,
                y = 16)) +
  
  # add in median absolute errors
  geom_text(data = filter(mae, timepoint == 5), 
            aes(label = paste0('MAE = ', mae), 
                x = 6.75, 
                y = 14)) +
  
  # formatting
  labs(title = 'Correlation and median absolute error (MAE) between epigenetic and chronological age ',
       x = 'Chronological Age',
       y = 'Epigenetic Age') +
  
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        plot.caption = element_text(hjust = 0),
        legend.position = 'bottom') +
  
  labs(caption = paste('Figure. Performance of epigenetic clocks from samples at follow-up.',
                       '\n',
                       'Dashed line is reference for perfect linear relationship',
                       '\n',
                       'r = Pearson correlation coefficient')
  )
followup

# output
################################################################################
# save plot as PNG
ggsave('figures/corr-mae-scatterplot-baseline.png', plot = baseline)
ggsave('figures/corr-mae-scatterplot-followup.png', plot = followup)