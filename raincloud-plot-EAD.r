### Author: Selene Banuelos adapted from code by Amy Inkster
### Date: 3/26/2026
### Description: Visualize distribution of epigenetic age deviation residuals 
### at baseline and follow-up as well as individual trajectories, colored by 
### PEARLS status (no/high)

# setup
library(tidyverse)
library(gghalves)
library(patchwork)

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
# function that generates EAD raincloud plot for given clock and tissue type
raincloud <- function(df, # data of interest
                      tissue, # (string) tissue type
                      clock # clock used to generate DNAm age for EAD
                      ){
  
  df %>%
  filter(Tissue == tissue) %>% # plot data from specified tissue type only
  ggplot(aes(x = Timepoint, # show EAD distribution at each timepoint
             y = !!sym(clock) # plot EAD from specified clock
             )
         ) +

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

  # formatting
  labs(
    y = 'EAD (years)',
    x = "Timepoint",
    color = "PEARLS"
    ) +
    
  scale_y_continuous(breaks = c(-2, 0, 2, 4),
                     limits = c(-2.5, 4.5)
                     ) +
    
  theme_classic() +
  theme(legend.position = 'bottom',
        # uncomment below for poster formatting
        legend.title = element_blank(),
        text = element_text(size = 24),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.subtitle = element_blank()
        )
  
}

# PedBE clock in blood samples
pbe_blood <- raincloud(clean_data, 'Blood', 'PedBEResid') +
  labs(title = 'Blood epigenetic age deviation (EAD) from chronological age across time',
       subtitle = 'PedBE')

# Skin & Blood clock in blood samples
sb_blood <- raincloud(clean_data, 'Blood', 'Horvath2Resid') +
  labs(subtitle = 'Skin & Blood') 

# add blood plots together
blood_plots <- pbe_blood + sb_blood

# PedBE clock in buccal samples
pbe_buccal <- raincloud(clean_data, 'Buccal', 'PedBEResid') +
  labs(title = "Buccal epigenetic age deviation (EAD) from chronological age across time",
       subtitle = 'PedBE')

# Skin & Blood clock in buccal samples
sb_buccal <- raincloud(clean_data, 'Buccal', 'Horvath2Resid') +
  labs(subtitle = 'Skin & Blood')

# add buccal plots together
buccal_plots <- pbe_buccal + sb_buccal

# output 
################################################################################
# save tissue-combined plots as PNG in figures folder
ggsave('figures/raincloud_ead_blood.png',
       plot = blood_plots,
       width = 7,
       height = 5,
       units = 'in')

ggsave('figures/raincloud_ead_buccal.png',
       plot = buccal_plots,
       width = 7,
       height = 5,
       units = 'in')

# save all plots individually
# ggsave('figures/raincloud-ead-horvath2-blood.png', 
#        plot = sb_blood,
#        width = 5.97,
#        height = 4.5,
#        units = 'in')
# 
# ggsave('figures/raincloud-ead-horvath2-buccal.png', 
#        plot = sb_buccal,
#        width = 5.97,
#        height = 4.5,
#        units = 'in')
# 
# ggsave('figures/raincloud-ead-pedbe-blood.png', 
#        plot = pbe_blood,
#        width = 5.97,
#        height = 4.5,
#        units = 'in')
# 
# ggsave('figures/raincloud-ead-pedbe-buccal.png', 
#        plot = pbe_buccal,
#        width = 5.97,
#        height = 4.5,
#        units = 'in')