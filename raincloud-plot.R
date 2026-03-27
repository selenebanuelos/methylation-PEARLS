### Author: Selene Banuelos adapted from code from Amy Inkster
### Date: 3/26/2026
### Description: Visualize distribution of epigenetic age deviation residuals 
### at baseline and follow-up as well as individual trajectories, colored by 
### PEARLS status (no/high)

# setup
library(tidyverse)
library(gghalves)

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
clean_data <- ead %>%
  # keep only EAD residuals
  select(SampleID, Horvath2Resid, PedBEResid) %>%
  # join residuals with sample information and PEARLS score
  full_join(sample_info, by = 'SampleID') %>%
  # create no PEARLS/high PEARLS groups from PEARLS score
  mutate(pearls = case_when(aces_baseline == 0 ~ 'no',
                            aces_baseline >= 5 ~ 'high')) %>%
  # change T2 to baseline and T5 to follow-up
  mutate(Timepoint = case_when(Timepoint == 'T2' ~ 'Baseline',
                               Timepoint == 'T5' ~ 'Follow-up'))

# data visualization
################################################################################
# Skin & Blood clock in blood samples
plot <- clean_data %>%
  filter(Tissue == 'Blood') %>%
  ggplot(aes(x = Timepoint, y = Horvath2Resid)) +
  
  # boxplot
  geom_boxplot(color = 'darkgrey', outlier.shape = NA, width = 0.1) +
  
  # paired lines that connect same-participant dots between time A and time B
  geom_line(aes(group = subjectid, color = pearls), alpha = 0.7) +
  
  # 1 - half violin density shading (baseline - left)
  gghalves::geom_half_violin(
    data = ~ filter(.x, Timepoint == 'Baseline'),
    fill = 'darkgrey',
    color = NA,
    side = "l",
    alpha = 0.6,
    nudge = 0.1,
    adjust = 0.5,
    trim = TRUE
  ) +
  
  # 2 - half violin, need 2nd call of this to do RHS of plot (follow-up - right)
  gghalves::geom_half_violin(
    data = ~ filter(.x, Timepoint == "Follow-up"),
    fill = 'darkgrey',
    color = NA,
    side = "r",
    alpha = 0.6,
    nudge = 0.1,
    adjust = 0.5,
    trim = TRUE
  ) +
  
  # individual participant data points, could add a geom_jitter to add X noise
  geom_point(aes(color = pearls), alpha = 0.6) +
  
  labs(
    title = "Epigenetic age deviation (EAD) from chronological age across time",
    subtitle = 'Skin & Blood clock in blood samples',
    y = "EAD (years)",
    x = "Timepoint",
    color = "Timepoint",
    fill = "Timepoint"
  ) +
  
  theme_classic() +
  theme(legend.position = "none")

plot