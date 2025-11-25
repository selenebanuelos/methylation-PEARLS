### Author: Selene Banuelos
### Date: 11/20/2025
### Description: Identify which samples did not pass QC

# setup
library(dplyr)
library(tidyr)

# import processed sample info
blood_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')
buccal_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# import master list of samples sent for processing
master <- read.csv('data-processed/pearls_acesmatchingbysexage.csv') %>%
  # keep only subject and specimen ids
  select(c(subjectid, T2_specimenid, T5_specimenid)) %>%
  # reshape data wide to long
  pivot_longer(
    cols = contains('specimenid'),
    names_to = 'Timepoint',
    values_to = 'specimenid',
    names_pattern = '(.+)_specimenid'
  ) 

# identify which blood samples didn't pass QC
all_blood <- select(blood_info, 
                          c(subjectid, specimenid, Timepoint, Tissue, SampleName)
                          ) %>%
  mutate(specimenid = as.numeric(as.character(specimenid))) %>%
  full_join(.,
            master,
            by = c('subjectid', 'specimenid', 'Timepoint')
  )

# identify which buccal samples didn't pass QC  
all_buccal <- select(buccal_info, 
                    c(subjectid, specimenid, Timepoint, Tissue, SampleName)
                    ) %>%
  mutate(specimenid = as.numeric(as.character(specimenid))) %>%
  full_join(.,
            master,
            by = c('subjectid', 'specimenid', 'Timepoint')
  )

# combine total dataset of samples that did and didn't pass QC after processing
all_samples <- rbind(all_blood, all_buccal)