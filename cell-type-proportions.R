### Author: Selene Banuelos
### Date: 1/21/2025
### Description: Visualize cell type proportions for each timepoint and tissue

# setup
library(dplyr)
library(tidyr)
library(ggplot2)
options(scipen = 999) # disable scientific notation

# import data
################################################################################
# cell type proportions estimated by Kobor lab
blood_cells <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds') %>%
  # only keep cell type proportion columns
  select(specimenid,Bas:Treg) %>%
  # calculate 'other' cell type proportion 
  mutate(other = 1 - rowSums(select(., !specimenid)))

buccal_cells <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds') %>%
  select(specimenid, Epi:Eosino) %>%
  mutate(other = 1 - rowSums(select(., !specimenid)))

# sample information
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv') %>%
  # transform to long format to get separate row for each timepoint
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
                      names_pattern = 'T(.)_(.*)', # remove 'T' from timepoint
                      names_to = c('timepoint', '.value')
                      ) %>%
  # create categorical ACEs variable (no = 0, high >= 5)
  mutate(aces = case_when(aces_baseline == 0 ~ 'no',
                          aces_baseline >= 5 ~ 'high')
  )

# data wrangling
################################################################################
# combine cell type proportions and sample timepoint/ACEs information 
blood <- sample %>%
  select(specimenid, timepoint, aces) %>%
  # mutate specimenid to factor for joining downstream
  mutate(specimenid = as.factor(specimenid)) %>%
  # join sample information with cell type proportions
  right_join(blood_cells,
             by = 'specimenid') %>%
  # make data longer for plotting
  pivot_longer(
    cols = Bas:other,
    names_to = 'cell_type',
    values_to = 'proportion'
  )

buccal <- sample %>%
  select(specimenid, timepoint, aces) %>%
  mutate(specimenid = as.factor(specimenid)) %>%
  right_join(buccal_cells,
             by = 'specimenid') %>%
  pivot_longer(
    cols = Epi:other,
    names_to = 'cell_type',
    values_to = 'proportion'
  )

# data visualization
################################################################################
# bar plots grouped by ACEs status
# blood samples timepoint 2
blood %>%
  filter(timepoint == '2') %>%
  ggplot(aes(x = cell_type,
             y = proportion,
             fill = aces)) +
  geom_boxplot() + 
  ggtitle('Blood samples T2')

# buccal samples timepoint 2
buccal %>%
  filter(timepoint == '2') %>%
  ggplot(aes(x = cell_type,
             y = proportion,
             fill = aces)) +
  geom_boxplot() +
  ggtitle('Buccal samples T2')

# blood samples timepoint 5
blood %>%
  filter(timepoint == '5') %>%
  ggplot(aes(x = cell_type,
             y = proportion,
             fill = aces)) +
  geom_boxplot() +
  ggtitle('Blood samples T5')

# buccal samples timepoint 5
buccal %>%
  filter(timepoint == '5') %>%
  ggplot(aes(x = cell_type,
             y = proportion,
             fill = aces)) +
  geom_boxplot() +
  ggtitle('Buccal samples T5')
