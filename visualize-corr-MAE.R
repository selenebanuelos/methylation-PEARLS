### Author: Selene Banuelos
### Date: 1/22/2026
### Description: Create scatterplots to visualize Pearson correlations between 
### chronological age and DNAm age generated using various clocks
### Calculate median absolute error and add to plots
### Correlations and MAE caclulated for each timepoint separately

# setup
library(dplyr)
library(ggplot2)
library(Metrics) # calculate meadian absolute error

# import data
################################################################################
DNAm_age <- read.csv('data-processed/dnam-age-sample-info.csv')

# data wrangling
################################################################################
# DNAm biomarkers of interest (predicting chronological age)
biomarkers <- c('horvath_age', 
                'horvath2', 
                'hannum_age', 
                'epigenetic_age_zhang', 
                'c_age', 
                'ped_be',
                'pheno_age',
                #'pace',
                'pc_horvath1',
                'pc_horvath2',
                'pc_hannum',
                'pc_pheno_age'
                )

# vector that only contains clocks with PC versions
pc <- c('horvath_age', 
        'pc_horvath1',
        'horvath2', 
        'pc_horvath2',
        'hannum_age', 
        'pc_hannum',
        'pheno_age',
        'pc_pheno_age')

# clean up DNAm and chrono age data for plotting
ages_long <- DNAm_age %>%
  select(-c(age_baseline, sex, imp_method, specimenid)) %>%
  #select(pearls_id, biomarkers, age, tissue, aces_baseline, timepoint) %>%
  # make data longer for plotting downstream
  tidyr::pivot_longer(
    cols = biomarkers,
    names_to = 'clock',
    values_to = 'dnam_age'
  ) %>%
  # create categorical PEARLS variable
  mutate(Pearls = case_when(aces_baseline == 0 ~ 'No',
                            aces_baseline >= 5 ~ 'High'
  ),
  # # format clock names for plot
  # clock = case_when(clock == 'horvath_age' ~ 'Horvath',
  #                   clock == 'horvath2' ~ 'Skin & Blood',
  #                   clock == 'hannum_age' ~ 'Hannum',
  #                   clock == 'ped_be' ~ 'PedBE',
  #                   clock == 'epigenetic_age_zhang' ~ 'Zhang',
  #                   clock == 'c_age' ~ 'cAge'),
  # format tissue names for plot
  tissue = case_when(tissue == 'blood' ~ 'Blood',
                     tissue == 'buccal' ~ 'Buccal'),
  # indicate if PC version of clock
  pc_version = case_when(grepl('pc_', clock) ~ 'yes',
                 .default = 'no')
  )
  

# calculate DNAm age ~ chrono age correlations 
################################################################################
cor <- ages_long %>%
  # calculate correlations at each timepoint separately for each tissue/clock group
  group_by(timepoint, tissue, clock, pc_version) %>%
  summarize(corr = round(cor(age, dnam_age),
                         digits = 2
                         )
            ) %>%
  ungroup()

# calculate DNAm age - chrono age median absolute error (MAE)
################################################################################
mae <- ages_long %>%
  # calculate MAE within each timepoint separately for each tissue/clock group
  group_by(timepoint, tissue, clock, pc_version) %>%
  summarize(mae = round(mdae(actual = age, predicted = dnam_age),
                        digits = 2)
            ) %>%
  ungroup()

# data visualization
################################################################################
# barplots to visualize median absolute error
# create named vector for facet labels in bar plots
label_names <- c(`2` = 'Baseline',
                 `5` = 'Follow-up',
                 `Blood` = 'Blood',
                 `Buccal` = 'Buccal')

# compare PC versions of clocks with their corresponding non-PC version
mae_pc_nonpc <- mae %>%
  filter(clock %in% pc) %>%
  ggplot(aes(x = clock,
             y = mae, 
             fill = pc_version)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = mae), position = position_dodge(width = 0.9), vjust = -0.25) +
  # order bars in the same order as elements of 'pc' vector
  scale_x_discrete(limits = pc) +
  # facet plots
  facet_grid(timepoint ~ tissue, labeller = as_labeller(label_names)) +
  # formatting
  labs(title = 'Median absolute error between epigenetic and chronological age',
       x = 'Clock',
       y = 'Median absolute error') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')
# non-PC versions have lower MAE across all clocks

# compare all non-PC version clocks
mae_nonpc <- mae %>%
  filter(pc_version == 'no') %>%
  ggplot(aes(x = clock,
             y = mae, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = mae), position = position_dodge(width = 0.9), vjust = -0.25) +
  # facet plots
  facet_grid(timepoint ~ tissue, labeller = as_labeller(label_names)) +
  # formatting
  labs(title = 'Median absolute error between epigenetic and chronological age',
       x = 'Clock',
       y = 'Median absolute error') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

# barplots to visualize Pearson correlation
################################################################################
# compare PC versions of clocks with their corresponding non-PC version
cor_pc_nonpc <- cor %>%
  filter(clock %in% pc) %>%
  ggplot(aes(x = clock,
             y = corr, 
             fill = pc_version)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = corr), position = position_dodge(width = 0.9), vjust = -0.25) +
  # order bars in the same order as elements of 'pc' vector
  scale_x_discrete(limits = pc) +
  # facet plots
  facet_grid(timepoint ~ tissue, labeller = as_labeller(label_names)) +
  # formatting
  labs(title = 'Pearson correlation coefficient between epigenetic and chronological age',
       x = 'Clock',
       y = 'r') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')
# what do you see?

# compare all non-PC version clocks
cor_nonpc <- cor %>%
  filter(pc_version == 'no') %>%
  ggplot(aes(x = clock,
             y = corr, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = corr), position = position_dodge(width = 0.9), vjust = -0.25) +
  # facet plots
  facet_grid(timepoint ~ tissue, labeller = as_labeller(label_names)) +
  # formatting
  labs(title = 'Pearson correlation coefficient between epigenetic and chronological age',
       x = 'Clock',
       y = 'r') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

# barplots for MAE and correlation for 3 biomarkers selected
mae_3 <- mae %>%
  filter(pc_version == 'no',
         clock %in% c('horvath_age', 'horvath2', 'ped_be')) %>%
  ggplot(aes(x = clock,
             y = mae, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = mae), position = position_dodge(width = 0.9), vjust = -0.25) +
  # facet plots
  facet_grid(timepoint ~ tissue, labeller = as_labeller(label_names)) +
  # formatting
  labs(title = 'Median absolute error between epigenetic and chronological age',
       x = 'Clock',
       y = 'Median absolute error') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

cor_3 <- cor %>%
  filter(pc_version == 'no',
         clock %in% c('horvath_age', 'horvath2', 'ped_be')) %>%
  ggplot(aes(x = clock,
             y = corr, 
             fill = clock)) +
  # create bar plot
  geom_bar(stat = "Identity") +
  # show MAE above each bar
  geom_text(aes(label = corr), position = position_dodge(width = 0.9), vjust = -0.25) +
  # facet plots
  facet_grid(timepoint ~ tissue, labeller = as_labeller(label_names)) +
  # formatting
  labs(title = 'Pearson correlation coefficient between epigenetic and chronological age',
       x = 'Clock',
       y = 'r') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold', size = 12),
        strip.background = element_rect(fill = 'darkgrey'),
        legend.position = 'none')

# output
################################################################################
# save plot as PNG
# ggsave('figures/corr-mae-scatterplot-baseline.png', plot = baseline)
# ggsave('figures/corr-mae-scatterplot-followup.png', plot = followup)

# MAE plots
ggsave('figures/mae-pc-nonpc.png', 
       plot = mae_pc_nonpc,
       width = 14,
       height = 9,
       units = 'in')
ggsave('figures/mae-nonpc.png', 
       plot = mae_nonpc,
       width = 16,
       height = 9,
       units = 'in')

# correlation plots
ggsave('figures/cor-pc-nonpc.png', 
       plot = cor_pc_nonpc,
       width = 14,
       height = 9,
       units = 'in')
ggsave('figures/cor-nonpc.png', 
       plot = cor_nonpc,
       width = 16,
       height = 9,
       units = 'in')