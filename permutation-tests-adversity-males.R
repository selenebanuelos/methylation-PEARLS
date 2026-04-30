### Author: Selene Banuelos adapted code by Mike Marin
### Date: 4/29/2026
### Description: Permutation test comparing mean EAD between adversity groups in
### male participants only

# setup
library(dplyr)
library(tidyr)
library(ggplot2)
options(scipen = 999)

# import data 
################################################################################
# epigenetic age deviation (EAD) residuals
ead <- read.csv('data-processed/methylation-predictions.csv')

# import processed sample info
blood_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')
buccal_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# participant ages at each visit
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
                            aces_baseline >= 5 ~ 'high')) %>%
  # this set of test are being done only in male participants
  filter()

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

# create function that performs permutation test
################################################################################
permutation_test <- function(df, # (dataframe) data
                             response, # (string) response variable to shuffle
                             p_n = 100000 # number of permutation samples to take
){
  
  # generate permutation samples ###############################################
  # the number of observations to sample (same size as original sample)
  n <- nrow(df)
  
  # number of permutation samples to take
  p <- 100000 
  
  # vector of the variable we're shuffling/sampling from (response variable)
  # could shuffle the labels (PEARLS group) instead
  var <- pull(df, response)
  
  # initialize a matrix (with all zeros) to store the permutation data
  perm_samples <- matrix(0, nrow = n, ncol = p_n) # each col is a permutation sample 
  
  # generate p_n permutation samples
  for (i in 1:p_n){
    
    # use sample() to take samples without replacement from original dataset
    perm_samples[,i] <- sample(var, size = n, replace = FALSE)
    
  }
  
  # calculate test stats for each permutation sample ###########################
  # initialize vector to store all of the test-stats
  perm_test_stats <- rep(0, p_n) # vector same length as number of permutation samples
  
  # loop through and calculate the test statistics
  for (i in 1:p_n){
    
    # calculate test-stat: difference in mean EAD between PEARLS groups
    perm_test_stats[i] <- mean(perm_samples[df$pearls == 'high', i]) -
      mean(perm_samples[df$pearls == 'no', i])
    # mean EAD in high PEARLS participants - mead EAD in no PEARLS participants
  }
  
  # calculate permutation p-value ##############################################
  # what is the probability of getting the observed test statistic or more 
  # extreme value if the null hypothesis of no difference is true?
  
  # get observed mean EAD in no PEARLS participants 
  mean_no <- filter(df, pearls == 'no') %>%
    pull(response) %>%
    mean()
  
  # get observed mean EAD in high PEARLS participants
  mean_high <- filter(df, pearls == 'high') %>%
    pull(response) %>%
    mean()
  
  # calculate observed test statistic 
  obs_test_stat <- mean_high - mean_no
  
  # p-value = # of perm test-stats >= observed test stat/ total # perm test-stats
  p_value <- sum(perm_test_stats >= obs_test_stat) / p_n
  
  return(sapply(c(obs_test_stat, p_value), round, digits = 2))
  
}

# conduct permutation tests in blood samples
################################################################################
# Skin & Blood clock, baseline
set.seed(123) # set seed for reproducibility
bl_sb_baseline <- filter(combined, 
                         Tissue == 'Blood',
                         Timepoint == 'T2') %>%
  permutation_test(., 'Horvath2Resid')

# Skin & Blood clock, follow-up
set.seed(123)
bl_sb_followup <- filter(combined, 
                         Tissue == 'Blood',
                         Timepoint == 'T5') %>%
  permutation_test(., 'Horvath2Resid')

# Skin & Blood clock, followup-baseline
set.seed(123)
bl_sb_diff <- filter(ead_diff,
                     Tissue == 'Blood',
                     !is.na(horvath2_diff)) %>%
  permutation_test(., 'horvath2_diff')

# SKin & Blood clock, EAD trajectory
set.seed(123)
bl_sb_traj <- filter(ead_diff,
                     Tissue == 'Blood',
                     !is.na(horvath2_trajectory)) %>%
  permutation_test(., 'horvath2_trajectory')

# PedBE clock, baseline
set.seed(123) 
bl_pbe_baseline <- filter(combined, 
                          Tissue == 'Blood',
                          Timepoint == 'T2') %>%
  permutation_test(., 'PedBEResid')

# PedBE clock, follow-up
set.seed(123) 
bl_pbe_followup <- filter(combined, 
                          Tissue == 'Blood',
                          Timepoint == 'T5') %>%
  permutation_test(., 'PedBEResid')

# PedBE clock, followup-baseline
set.seed(123)
bl_pbe_diff <- filter(ead_diff,
                      Tissue == 'Blood',
                      !is.na(pedbe_diff)) %>%
  permutation_test(., 'pedbe_diff')

# PedBE clock, EAD trajectory
set.seed(123)
bl_pbe_traj <- filter(ead_diff,
                      Tissue == 'Blood',
                      !is.na(pedbe_trajectory)) %>%
  permutation_test(., 'pedbe_trajectory')

# conduct permutation tests in buccal samples
################################################################################
# Skin & Blood clock, baseline
set.seed(123) # set seed for reproducibility
bu_sb_baseline <- filter(combined, 
                         Tissue == 'Buccal',
                         Timepoint == 'T2') %>%
  permutation_test(., 'Horvath2Resid')

# Skin & Blood clock, follow-up
set.seed(123)
bu_sb_followup <- filter(combined, 
                         Tissue == 'Buccal',
                         Timepoint == 'T5') %>%
  permutation_test(., 'Horvath2Resid')

# Skin & Blood clock, followup-baseline
set.seed(123)
bu_sb_diff <- filter(ead_diff,
                     Tissue == 'Buccal',
                     !is.na(horvath2_diff)) %>%
  permutation_test(., 'horvath2_diff')

# Skin & Blood clock, EAD trajectory
set.seed(123)
bu_sb_traj <- filter(ead_diff,
                     Tissue == 'Buccal',
                     !is.na(horvath2_trajectory)) %>%
  permutation_test(., 'horvath2_trajectory')

# PedBE clock, baseline
set.seed(123) 
bu_pbe_baseline <- filter(combined, 
                          Tissue == 'Buccal',
                          Timepoint == 'T2') %>%
  permutation_test(., 'PedBEResid')

# PedBE clock, follow-up
set.seed(123) 
bu_pbe_followup <- filter(combined, 
                          Tissue == 'Buccal',
                          Timepoint == 'T5') %>%
  permutation_test(., 'PedBEResid')

# PedBE clock, followup-baseline
set.seed(123)
bu_pbe_diff <- filter(ead_diff,
                      Tissue == 'Buccal',
                      !is.na(pedbe_diff)) %>%
  permutation_test(., 'pedbe_diff')

# PedBE clock, EAD trajectory
set.seed(123)
bu_pbe_traj <- filter(ead_diff,
                      Tissue == 'Buccal',
                      !is.na(pedbe_trajectory)) %>%
  permutation_test(., 'pedbe_trajectory')