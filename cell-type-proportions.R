### Author: Selene Banuelos
### Date: 1/21/2025
### Description: Visualize cell type proportions for each timepoint and tissue

# setup
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
options(scipen = 999) # disable scientific notation

# import data
################################################################################
# cell type proportions estimated by Kobor lab
blood_cells <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds') %>%
  # only keep cell type proportion columns
  select(specimenid,Bas:Treg)

buccal_cells <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds') %>%
  select(specimenid, Epi:Eosino) %>%
  # calculate proportion of immune cells
  mutate(IC = rowSums(select(., B, NK, CD4T, CD8T, Mono, Neutro, Eosino))
         )

# sample information
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv') %>%
  # transform to long format to get separate row for each timepoint
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
                      names_pattern = 'T(.)_(.*)', # remove 'T' from timepoint
                      names_to = c('timepoint', '.value')
                      ) %>%
  # create categorical ACEs variable (no = 0, high >= 5)
  mutate(aces = case_when(aces_baseline == 0 ~ 'no',
                          aces_baseline >= 5 ~ 'high'),
         sex = as.character(sex)
         )

# data wrangling
################################################################################
# combine cell type proportions and sample timepoint/ACEs information 
blood <- sample %>%
  select(specimenid, timepoint, aces, sex) %>%
  # mutate specimenid to factor for joining downstream
  mutate(specimenid = as.factor(specimenid)) %>%
  # join sample information with cell type proportions
  right_join(blood_cells,
             by = 'specimenid') %>%
  # make data longer for plotting
  pivot_longer(
    cols = Bas:Treg,
    names_to = 'cell_type',
    values_to = 'proportion'
  ) %>%
  # calculate average proportion of each cell type for each adversity group
  group_by(aces, cell_type, timepoint) %>%
  mutate(mean_prop = mean(proportion)) %>%
  ungroup()

# buccal samples including all cell types
buccal <- sample %>%
  select(specimenid, timepoint, aces, sex) %>%
  mutate(specimenid = as.factor(specimenid)) %>%
  right_join(buccal_cells,
             by = 'specimenid') %>%
 pivot_longer(
    cols = Epi:Eosino,
    names_to = 'cell_type',
    values_to = 'proportion'
  ) %>%
  group_by(aces, cell_type, timepoint) %>%
  mutate(mean_prop = mean(proportion)) %>%
  ungroup()

# buccal samples separated into 3 types: epithelial, fibroblasts, immune cells
buccal_3 <- sample %>%
  select(specimenid, timepoint, aces, sex) %>%
  mutate(specimenid = as.factor(specimenid)) %>%
  right_join(buccal_cells,
             by = 'specimenid') %>%
  rowwise(specimenid) %>%
  # remove immune cell subtypes
  select(-c(B, NK, CD4T, CD8T, Mono, Neutro, Eosino)) %>%
  pivot_longer(
    cols = c(Epi, Fib, IC),
    names_to = 'cell_type',
    values_to = 'proportion'
  ) %>%
  group_by(aces, cell_type, timepoint) %>%
  mutate(mean_prop = mean(proportion)) %>%
  ungroup()

# buccal samples separated into 7 immune cell subtypes
buccal_7 <- sample %>%
  select(specimenid, timepoint, aces, sex) %>%
  mutate(specimenid = as.factor(specimenid)) %>%
  right_join(buccal_cells,
             by = 'specimenid') %>%
  rowwise(specimenid) %>%
  # remove immune cell subtypes
  select(-c(Epi, Fib, IC)) %>%
  pivot_longer(
    cols = B:Eosino,
    names_to = 'cell_type',
    values_to = 'proportion'
  ) %>%
  group_by(aces, cell_type, timepoint) %>%
  mutate(mean_prop = mean(proportion)) %>%
  ungroup()

# data visualization
################################################################################
### blood samples
# side by side bar plots
bar_bl <- blood %>%
  ggplot(aes(x = cell_type, 
             y = proportion, 
             fill = cell_type, 
             color = cell_type, 
             alpha = aces # vary opaqueness by adversity group
             )) +
  geom_bar(stat = "summary_bin", fun = mean, position = "dodge") +
  scale_alpha_manual(values = c("no" = 0.25, "high" = 1)) +
  facet_grid(. ~timepoint)

# box plots stratified by adversity status and shown for each timepoint
box_bl <- blood %>%
  ggplot(aes(x = cell_type, y = proportion, fill = aces)) +
  geom_boxplot() +
  facet_grid(. ~timepoint)

# spaghetti plot
spaghet_bl <- blood %>%
  ggplot(aes(x = timepoint,y = mean_prop, group = cell_type, color = cell_type)) +
  geom_point() + 
  geom_line() +
  facet_grid(. ~aces)

bl <- bar_bl + box_bl + spaghet_bl + plot_layout(ncol = 1)

### buccal samples, including all cell types
# side by side bar plots
bar_buc <- buccal %>%
  ggplot(aes(x = cell_type, 
             y = proportion, 
             fill = cell_type, 
             color = cell_type, 
             alpha = aces # vary opaqueness by adversity group
  )) +
  geom_bar(stat = "summary_bin", fun = mean, position = "dodge") +
  scale_alpha_manual(values = c("no" = 0.25, "high" = 1)) +
  facet_grid(. ~timepoint)

# box plots stratified by adversity status and shown for each timepoint
box_buc <- buccal %>%
  ggplot(aes(x = cell_type, y = proportion, fill = aces)) +
  geom_boxplot() +
  facet_grid(. ~timepoint)

# spaghetti plot
spaghet_buc <- buccal %>%
  ggplot(aes(x = timepoint,y = mean_prop, group = cell_type, color = cell_type)) +
  geom_point() + 
  geom_line() +
  facet_grid(. ~aces)

buc <- bar_buc + box_buc + spaghet_buc + plot_layout(ncol = 1)

### buccal samples: epi/fib/IC
# side by side bar plots
bar_buc_3 <- buccal_3 %>%
  ggplot(aes(x = cell_type, 
             y = proportion, 
             fill = cell_type, 
             color = cell_type, 
             alpha = aces # vary opaqueness by adversity group
  )) +
  geom_bar(stat = "summary_bin", fun = mean, position = "dodge") +
  scale_alpha_manual(values = c("no" = 0.25, "high" = 1)) +
  facet_grid(. ~timepoint)

# box plots stratified by adversity status and shown for each timepoint
box_buc_3 <- buccal_3 %>%
  ggplot(aes(x = cell_type, y = proportion, fill = aces)) +
  geom_boxplot() +
  facet_grid(. ~timepoint)

# spaghetti plot
spaghet_buc_3 <- buccal_3 %>%
  ggplot(aes(x = timepoint,y = mean_prop, group = cell_type, color = cell_type)) +
  geom_point() + 
  geom_line() +
  facet_grid(. ~aces)

buc_3 <- bar_buc_3 + box_buc_3 + spaghet_buc_3 + plot_layout(ncol = 1)

### buccal samples: IC subtypes
# side by side bar plots
bar_buc_7 <- buccal_7 %>%
  ggplot(aes(x = cell_type, 
             y = proportion, 
             fill = cell_type, 
             color = cell_type, 
             alpha = aces # vary opaqueness by adversity group
  )) +
  geom_bar(stat = "summary_bin", fun = mean, position = "dodge") +
  scale_alpha_manual(values = c("no" = 0.25, "high" = 1)) +
  facet_grid(. ~timepoint)

# box plots stratified by adversity status and shown for each timepoint
box_buc_7 <- buccal_7 %>%
  ggplot(aes(x = cell_type, y = proportion, fill = aces)) +
  geom_boxplot() +
  facet_grid(. ~timepoint)

# spaghetti plot
spaghet_buc_7 <- buccal_7 %>%
  ggplot(aes(x = timepoint,y = mean_prop, group = cell_type, color = cell_type)) +
  geom_point() + 
  geom_line() +
  facet_grid(. ~aces)

buc_7 <- bar_buc_7 + box_buc_7 + spaghet_buc_7 + plot_layout(ncol = 1)

# output
################################################################################
ggsave('figures/cell-props-blood.png', bl, width = 11, height = 11, units = 'in')
ggsave('figures/cell-props-buccal.png', buc, width = 11, height = 11, units = 'in')
ggsave('figures/cell-props-buccal-3.png', buc_3, width = 11, height = 11, units = 'in')
ggsave('figures/cell-props-buccal-7.png', buc_7, width = 11, height = 11, units = 'in')