### Author: Selene Banuelos
### Date: 1/21/2025
### Description: Visualize cell type proportions for each timepoint and tissue

# setup
library(dplyr)
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
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# data wrangling
################################################################################
# combine cell type proportions and sample timepoint/tissue information 

# data visualization
################################################################################
# bar plots
# blood samples timepoint 2

# buccal samples timepoint 2

# blood samples timepoint 5

# buccal samples timepoint 5
