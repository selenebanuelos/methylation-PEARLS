### Author: Selene Banuelos
### Date: 1/22/2026
### Description: Create scatterplots to visualize correlations between DNAm age
### and chronological age for predictions generated using PedBE & Horvath2 clocks
### Calculate median absolute error and add to plots

# setup
library(dplyr)
library(ggplot2)
library(rmcorr) # calculate correlation for repeated measures
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
# buccal samples, PedBE clock
bu_pbe <- ages_long %>%
  filter(tissue == 'Buccal', 
         clock == 'PedBE'
         ) %>%
  rmcorr(participant = pearls_id,
                measure1 = age,
                measure2 = dnam_age,
                dataset = .
                )

# buccal samples, Skin & Blood clock
bu_sb <- ages_long %>%
  filter(tissue == 'Buccal', 
         clock == 'Skin & Blood'
  ) %>%
  rmcorr(participant = pearls_id, 
         measure1 = age, 
         measure2 = dnam_age, 
         dataset = .
  )

# blood samples, PedBE clock
bl_pbe <- ages_long %>%
  filter(tissue == 'Blood', 
         clock == 'PedBE'
  ) %>%
  rmcorr(participant = pearls_id, 
         measure1 = age, 
         measure2 = dnam_age, 
         dataset = .
  )

# blood samples, Skin & Blood clock
bl_sb <- ages_long %>%
  filter(tissue == 'Blood', 
         clock == 'Skin & Blood'
  ) %>%
  rmcorr(participant = pearls_id, 
         measure1 = age, 
         measure2 = dnam_age, 
         dataset = .
  )

# combine correlations for each tissue/clock group
corrs <- data.frame(
  tissue = c('Buccal', 'Buccal', 'Blood', 'Blood'),
  clock = c('PedBE', 'Skin & Blood', 'PedBE', 'Skin & Blood'),
  corr = c(bu_pbe$r, bu_sb$r, bl_pbe$r, bl_sb$r)
  ) %>%
  mutate(corr = round(corr, digits = 2)) %>%
  # join correlation coefficients to rest of the data
  right_join(ages_long, by = c('tissue', 'clock'))

# calculate DNAm age - chrono age median absolute error 
################################################################################
mae <- ages_long %>%
  # calculate MAE within each clock and tissue group
  group_by(tissue, clock) %>%
  summarize(mae = round(mdae(actual = age, predicted = dnam_age),
                        digits = 2
                        )
            ) %>%
  ungroup()

# data visualization
################################################################################
# create scatter plots with regression line to display correlation between
# epigenetic age and chronological age across both timepoints for each tissue
plot <- corrs %>%
  # only plot data from blood samples
  ggplot(aes(x = age,
             y = dnam_age
             )
         ) +
  
  # scatterplot showing relationship between DNAm age and chrono age
  geom_point(shape = 16, aes(color = Pearls)) + # color code points by ACEs status
  
  # add in repeated measures correlation line
  #geom_abline(aes(slope = corr)) +
  
  # stratify plot by epigenetic clock
  facet_grid(clock ~ tissue) +
  
  # add in reference line showing 100% correlation
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'darkgrey') +
  
  # add in correlation coefficients (r) and median absolute errors
  geom_text(data = corrs, aes(label = paste0('r_rm = ', corr), 
                              x = 2.4, 
                              y = 16) 
            ) +
  geom_text(data = mae, aes(label = paste0('MAE = ', mae)), x = 3, y = 14) +
  
  # formatting
  labs(title = 'Correlation and median absolute error (MAE) between epigenetic and chronological age ',
       x = 'Chronological Age',
       y = 'Epigenetic Age') +
  
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        plot.caption = element_text(hjust = 0),
        legend.position = 'bottom') +
  
  labs(caption = paste(' Dashed line is reference for perfect linear relationship',
                       '\n',
                       'r_rm = repeated measures correlation coefficient'
                       )
  )
plot

# output
################################################################################
# save plot as PNG
ggsave('figures/scatterplot-corr-mae.png',
       plot = plot)