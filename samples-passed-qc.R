### Author: Selene Banuelos
### Date: 11/20/2025
### Description: Identify which samples did not pass QC

# setup
library(dplyr)
library(tidyr)

# import data
################################################################################
# processed DNAme sample info
blood_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')
buccal_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# import master list of samples sent for processing
master <- read.csv('data-processed/pearls-acesmatchingbysexage.csv') %>%
  # keep only subject and specimen ids
  select(c(subjectid, T2_specimenid, T5_specimenid)) %>%
  # reshape data wide to long
  pivot_longer(
    cols = contains('specimenid'),
    names_to = 'Timepoint',
    values_to = 'specimenid',
    names_pattern = '(.+)_specimenid'
  ) 

# data wrangling
################################################################################
# identify which blood samples didn't pass QC
all_blood <- blood_info %>% # all DNAme samples that passed QC
  select(subjectid, specimenid, Timepoint) %>%
  # change specimenid to numeric for joining downstream
  mutate(specimenid = as.numeric(as.character(specimenid)),
         # create var that specifies sample passed QC (1 = yes)
         passed_qc = 1
         ) %>%
  # join with master list of samples sent in for processing
  full_join(.,
            master,
            by = c('subjectid', 'specimenid', 'Timepoint')
            ) %>%
  mutate(tissue = 'blood', # add tissue type var back in
         # identify samples that didnt pass qc with '0'
         passed_qc = case_when(is.na(passed_qc) ~ 0,
                               .default = passed_qc)
         )

# identify which buccal samples didn't pass QC  
all_buccal <- buccal_info %>% # all DNAme samples that passed QC
  select(subjectid, specimenid, Timepoint) %>%
  # change specimenid to numeric for joining downstream
  mutate(specimenid = as.numeric(as.character(specimenid)),
         # create var that specifies sample passed QC (1 = yes)
         passed_qc = 1
  ) %>%
  # join with master list of samples sent in for processing
  full_join(.,
            master,
            by = c('subjectid', 'specimenid', 'Timepoint')
  ) %>%
  mutate(tissue = 'buccal', # add tissue type var back in
         # identify samples that didnt pass qc with '0'
         passed_qc = case_when(is.na(passed_qc) ~ 0,
                               .default = passed_qc)
  )

# combine total dataset of samples that did and didn't pass QC after processing
all_samples <- rbind(all_blood, all_buccal)

# output 
################################################################################
# save csv with data on DNAme samples that did/didn't pass QC
write.csv(all_samples,
          file = 'data-processed/samples-passed-qc.csv',
          row.names = FALSE)