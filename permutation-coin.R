### Author: Selene Banuelos
### Date: 5/21/2026
### Description: Permutation test comparing mean EAD between adversity groups
### Testing entire sample and doing sex-stratified tests
### Doing exact permutations with coin package

# setup
library(dplyr)
library(tidyr)
library(coin)
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
# list of biomarkers to do permutation testing on
biomarkers <- list('HorvathAgeResid', 'Horvath2Resid', 'PedBEResid')

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
  select(SampleID, HorvathAgeResid, Horvath2Resid, PedBEResid) %>%
  # join residuals with sample information and PEARLS score
  full_join(sample_info, by = 'SampleID') %>%
  # create no PEARLS/high PEARLS groups from PEARLS score
  mutate(pearls = case_when(aces_baseline == 0 ~ 'no',
                            aces_baseline >= 5 ~ 'high'),
         pearls = as.factor(pearls)) 

# stratify dataset by biological sex
male <- filter(combined, sex_.0.F. == 1) # 1 = male
female <- filter(combined, sex_.0.F. == 0) # 0 = female

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

# define functions
################################################################################
# function to do exact permutation testing using coin::oneway_test
permutation_test <- function(outcome, # DNAm estimate name as string
                             tissue,
                             timepoint,
                             data
) {
  # filter data for specified tissue and timepoint
  df <- data %>%
    filter(Tissue == tissue,
           Timepoint == timepoint)
  
  # get outcome column as vector
  y <- df[[outcome]]
  
  # get group column as vector
  group <- df$pearls
  
  # exact permutation test
  result <- coin::oneway_test(formula = y ~ group, 
                              distribution = 'exact',
                              alternative = 'two.sided')
  
  # get observed mean outcome in no PEARLS participants 
  mean_no <- filter(df, pearls == 'no') %>%
    pull(outcome) %>%
    mean()
  
  # get observed mean outcome in high PEARLS participants
  mean_high <- filter(df, pearls == 'high') %>%
    pull(outcome) %>%
    mean()
  
  # calculate observed test statistic 
  obs_test_stat <- abs(mean_high - mean_no)

  
  # return data frame with results
  data.frame(tissue = tissue, # tissue type
             clock = outcome,
             timepoint = timepoint,
             test_stat = obs_test_stat,
             pvalue = pvalue(result), # get p-value from permutation test
             data = deparse(substitute(data)) # dataframe name
  )
}

# exact permutation testing
################################################################################
# do permutation tests for each biomarker in list
blood_2 <- map_df(biomarkers, permutation_test, 'Blood', 'T2', combined)
blood_5 <- map_df(biomarkers, permutation_test, 'Blood', 'T5', combined)

buccal_2 <- map_df(biomarkers, permutation_test, 'Buccal', 'T2', combined)
buccal_5 <- map_df(biomarkers, permutation_test, 'Buccal', 'T5', combined)

# sex-stratified testing 
m_blood_2 <- map_df(biomarkers, permutation_test, 'Blood', 'T2', male)
m_blood_5 <- map_df(biomarkers, permutation_test, 'Blood', 'T5', male)

m_buccal_2 <- map_df(biomarkers, permutation_test, 'Buccal', 'T2', male)
m_buccal_5 <- map_df(biomarkers, permutation_test, 'Buccal', 'T5', male)

f_blood_2 <- map_df(biomarkers, permutation_test, 'Blood', 'T2', female)
f_blood_5 <- map_df(biomarkers, permutation_test, 'Blood', 'T5', female)

f_buccal_2 <- map_df(biomarkers, permutation_test, 'Buccal', 'T2', female)
f_buccal_5 <- map_df(biomarkers, permutation_test, 'Buccal', 'T5', female)
