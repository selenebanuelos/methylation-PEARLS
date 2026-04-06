### Author: Selene Banuelos
### Date: 4/5/2025
### Description: Permutation test 

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

# data wrangling 
################################################################################
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

# generate permutation test statistics
################################################################################
# set seed for reproducibility
set.seed(123) 

# the number of observations to sample (same size as original sample)
n <- nrow(combined)

# number of permutation samples to take
p <- 100000 

# the variable we're shuffling (response variable here)
# could shuffle the labels (PEARLS group) instead
var <- combined$Horvath2Resid

# initialize a matrix to store the permutation data
perm_samples <- matrix(0, nrow = n, ncol = p) # each col is a permutation sample 

# generate p permutation samples
for (i in 1:p){
  
  # use sample() to take samples without replacement from original dataset
  perm_samples[,i] <- sample(var, size = n, replace = FALSE)
  
}

# use a loop to calculate test statistics for each sample (diff in means)
# initialize vectors to store all of the test-stats
test_stat <- rep(0, p) # vector same length as number of permutation samples

# loop through and calculate the test statistics
for (i in 1:p){
  
  # calculate test-stat: difference in mean EAD between PEARLS groups
  test_stat[i] <- mean(perm_samples[combined$pearls == 'no', i]) -
    mean(perm_samples[combined$pearls == 'high', i])
  # mean EAD in no PEARLS participants - mead EAD in high PEARLS participants
  
}

# calculate permutation p-value
################################################################################

# what is the probability of getting the observed test statistic (diff in mean EAD)
# or more extreme if the null hypothesis of no difference is true?

# get observed test statistic 
obs_test_stat <- mean(combined[combined$pearls == 'no', 'Horvath2Resid']) -
  mean(combined[combined$pearls == 'high', 'Horvath2Resid'])

# p-value = # of perm test-stats >= observed test stat/ total # perm test-stats
p_value <- sum(test_stat >= obs_test_stat) / length(test_stat)

# visualize distribution of permutation test statistics
################################################################################
perm_test_density <- with(density(test_stat), data.frame(x,y))

ggplot(perm_test_density, aes(x = x, y = y)) +
  geom_line() +
  geom_area(aes(x = ifelse(x > p_value, x, 0)), fill = 'red')


  geom_vline(xintercept = obs_test_stat, color = 'red', linewidth = 1)