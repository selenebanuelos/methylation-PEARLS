### Author: Selene Banuelos
### Date: 5/21/2026
### Description: Permutation test comparing mean EAD between adversity groups
### Testing entire sample and doing sex-stratified tests
### Doing exact permutations with coin package

# setup
library(dplyr)
library(tidyr)
library(coin)
library(ggplot2)
options(scipen = 999)

# import data 
################################################################################
# epigenetic age deviation (EAD) residuals
ead <- read.csv('data-processed/methylation-predictions.csv')

# import processed sample info
blood_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')
buccal_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# participant ages at each visit (to measure time between visits)
ages <- read.csv('data-raw/pearls_data_LauraDiaz_2025_11_20.csv') %>%
  select(pearls_id, collectionage_t2, collectionage_t5) %>%
  rename(subjectid = pearls_id)

# data wrangling 
################################################################################
# sample information variables
info_vars <- c('subjectid', 
               'Timepoint', 
               'Tissue', 
               'sex_.0.F.', 
               'aces_baseline', 
               'Sample_Name')

# combine sample information for blood and buccal samples
sample_info <- rbind(
  select(blood_info, info_vars),
  select(buccal_info, info_vars)
) %>%
  # rename sample ID for downstream joining
  rename(SampleID = Sample_Name)

# join aar data with sample info and clean up for analysis
combined <- ead %>%
  # keep only EAD residuals
  select(SampleID, Horvath2Resid, PedBEResid) %>%
  # join residuals with sample information and PEARLS score
  full_join(sample_info, by = 'SampleID') %>%
  # create no PEARLS/high PEARLS groups from PEARLS score
  mutate(pearls = case_when(aces_baseline == 0 ~ 'no',
                            aces_baseline >= 5 ~ 'high'),
         pearls = as.factor(pearls)) %>%
  # conducting permutation testing only in male participants 
  filter(sex_.0.F. == 1) # 1 = male

# calculate difference in EAD and EAD trajectory
ead_diff <- combined %>%
  pivot_wider(id_cols = c(subjectid, Tissue, pearls),
              names_from = Timepoint,
              values_from = c(Horvath2Resid, PedBEResid)
  ) %>%
  # add in subject age at each visit
  left_join(ages, by = 'subjectid') %>%
  # calculate difference in EAD
  mutate(horvath2_diff = (Horvath2Resid_T5 - Horvath2Resid_T2),
         pedbe_diff = (PedBEResid_T5 - PedBEResid_T2)
  ) %>%
  # calculate time between baseline and follow-up visits
  mutate(time_between = collectionage_t5 - collectionage_t2) %>%
  # calculate EAD trajectory (follow-up EAD - baseline EAD)/time between visits
  mutate(horvath2_trajectory = horvath2_diff/time_between,
         pedbe_trajectory = pedbe_diff/time_between)

blood_2 <- combined %>%
  filter(Tissue == 'Blood',
         Timepoint == 'T2')
blood_5 <- combined %>%
  filter(Tissue == 'Blood',
         Timepoint == 'T5')
buccal_2 <- combined %>%
  filter(Tissue == 'Buccal',
         Timepoint == 'T2')
buccal_5 <- combined %>%
  filter(Tissue == 'Buccal',
         Timepoint == 'T5')

all <- list(blood_2, blood_5, buccal_2, buccal_5)
# exact permutation testing
################################################################################
horvath2 <- lapply(all,
                   oneway_test(Horvath2Resid ~ pearls, 
                               data = ., 
                               distribution = approximate(nresample = 100000),
                               alternative = 'greater'
                               )
                   )

# approximate (what i did before)
# difference in mean outcome
oneway_test(Horvath2Resid ~ pearls, 
            data = test, 
            ytrafo = function(Horvath2Resid) as.numeric(Horvath2Resid[[1]]), 
            distribution = approximate(nresample = 100000),
            alternative = 'greater'
            )
# difference in median outcomes
oneway_test(Horvath2Resid ~ pearls, 
            data = test, 
            distribution = approximate(nresample = 100000),
            alternative = 'greater'
)

oneway_test(PedBEResid ~ pearls, data = test, distribution = 'exact', alternative = 'greater')