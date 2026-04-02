### Author: Selene Banuelos
### Date: 3/25/2026
### Description: Run t-tests to compare average epigenetic age deviation between
### no PEARLS and high PEARLS groups at baseline, follow-up, and difference
### between baseline and follow-up

# setup
library(dplyr)
library(tidyr)

# import data ##################################################################
# epigenetic age deviation (EAD) residuals
ead <- read.csv('data-processed/methylation-predictions.csv')

# import processed sample info
blood_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')
buccal_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# data wrangling ###############################################################
# sample information variables
info_vars <- c('subjectid', 'Timepoint', 'Tissue', 'aces_baseline', 'Sample_Name')

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
                            aces_baseline >= 5 ~ 'high'))

# calculate absolute difference in EAD
ead_diff <- combined %>%
  pivot_wider(id_cols = c(subjectid, Tissue, pearls),
              names_from = Timepoint,
              values_from = c(Horvath2Resid, PedBEResid)
              ) %>%
  # calculate absolute difference in EAD
  mutate(horvath2_diff = (Horvath2Resid_T5 - Horvath2Resid_T2),
         pedbe_diff = (PedBEResid_T5 - PedBEResid_T2)
  )

# define function ##############################################################
# function that extracts mean difference between groups and 95% CI
get_estimates <- function(results){ # t.test object
  
  # get mean outcome for high PEARLS group
  pearls_high <- results$estimate[[1]]
  
  # get mean outcome for no PEARLS group
  pearls_no <- results$estimate[[2]]
  
  # calculate difference in means
  diff <- round(pearls_high - pearls_no, digits = 2)
  
  # get 95% confidence intervals
  lower <- round(results$conf.int[1], digits = 2)
  upper <- round(results$conf.int[2], digits = 2)
  ci <- paste0('(', lower, ', ', upper, ')')
  
  # return difference in means and 95% CI
  return(c(diff, ci))
  
}

# t-tests in blood samples #####################################################
# compare mean Skin & Blood EAD between PEARLS groups at baseline
bl_sb_t2 <- t.test(Horvath2Resid ~ pearls, 
                   data = combined %>%
                     filter(Timepoint == 'T2',
                            Tissue == 'Blood')) %>%
  get_estimates()

# compare mean PedBE EAD between PEARLS groups at baseline
bl_pbe_t2 <- t.test(PedBEResid ~ pearls, 
                    data = combined %>%
                      filter(Timepoint == 'T2',
                             Tissue == 'Blood')) %>%
  get_estimates()

# compare mean Skin & Blood EAD between PEARLS groups at follow-up
bl_sb_t5 <- t.test(Horvath2Resid ~ pearls, 
                   data = combined %>%
                     filter(Timepoint == 'T5',
                            Tissue == 'Blood')) %>%
  get_estimates()

# compare mean PedBE EAD between PEARLS groups at follow-up
bl_pbe_t5 <- t.test(PedBEResid ~ pearls, 
                    data = combined %>%
                      filter(Timepoint == 'T5',
                             Tissue == 'Blood')) %>%
  get_estimates()

# compare mean Skin & Blood EAD diff across time between PEARLS groups
bl_sb_diff <- t.test(horvath2_diff ~ pearls, 
                     data = filter(ead_diff, Tissue == 'Blood')) %>%
  get_estimates()

# compare PedBE EAD diff across time between PEARLS groups
bl_pbe_diff <- t.test(pedbe_diff ~ pearls, 
                      data = filter(ead_diff, Tissue == 'Blood')) %>%
  get_estimates() 

# conduct t-tests in buccal samples #############################################
# compare mean Skin & Blood EAD between PEARLS groups at baseline
bu_sb_t2 <- t.test(Horvath2Resid ~ pearls, 
                   data = combined %>%
                     filter(Timepoint == 'T2',
                            Tissue == 'Buccal')) %>%
  get_estimates()

# compare mean PedBE EAD between PEARLS groups at baseline
bu_pbe_t2 <- t.test(PedBEResid ~ pearls, 
                    data = combined %>%
                      filter(Timepoint == 'T2',
                             Tissue == 'Buccal')) %>%
  get_estimates()

# compare mean Skin & Blood EAD between PEARLS groups at follow-up
bu_sb_t5 <- t.test(Horvath2Resid ~ pearls, 
                   data = combined %>% 
                     filter(Timepoint == 'T5',
                            Tissue == 'Buccal')) %>%
  get_estimates()

# compare mean PedBE EAD between PEARLS groups at follow-up
bu_pbe_t5 <- t.test(PedBEResid ~ pearls, 
                    data = combined %>%
                      filter(Timepoint == 'T5',
                             Tissue == 'Buccal')) %>%
  get_estimates()

# compare mean Skin & Blood EAD diff across time between PEARLS groups
bu_sb_diff <- t.test(horvath2_diff ~ pearls, 
                     data = filter(ead_diff, Tissue == 'Buccal')) %>%
  get_estimates()

# compare PedBE EAD diff across time between PEARLS groups
bu_pbe_diff <- t.test(pedbe_diff ~ pearls, 
       data = filter(ead_diff, Tissue == 'Buccal')) %>%
  get_estimates()