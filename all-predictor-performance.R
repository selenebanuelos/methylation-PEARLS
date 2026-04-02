### Author: Selene Banuelos
### Date: 12/15/2025
### Description: join methylation predictions with sample information

# setup
library(tidyverse)

# import data
################################################################################
# information on required CpGs per clock
clock_info <- read.csv('data-processed/summary_methscore_CpG.csv')

# generated methylation predictions
predictions <- read.csv('data-processed/methylation-predictions.csv')

# sample information
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# data wrangling 
################################################################################
# names of chronological age predictors
chrono_age <- c('horvath_age', 
                'horvath2', 
                'hannum_age', 
                'epigenetic_age_zhang', 
                'c_age', 
                'ped_be'
                )

# calculate % of missing predictive CpGs and filter for chronological age clocks
clocks <- clock_info %>%
  # format predictor name for joining downstream
  mutate(predictor = janitor::make_clean_names(predictor),
         # calculate % of missing CpGs
         perc_missing_cpg = 100*(nCpG_missing / nCpG_required)
  ) %>%
  filter(predictor %in% chrono_age)

# clean up DNAme sample info for joining with methylation predictions 
sample_long <- sample %>%
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
               names_pattern = '(.)_(.*)',
               names_to = c('timepoint', '.value')
               ) %>%
  mutate(specimenid = as.character(specimenid)) %>%
  dplyr::rename(pearls_id = subjectid)

# standardize format of clock names to match chrono_age vec and for joining
names(predictions) <- janitor::make_clean_names(names(predictions))

# join methylation predictions with sample information
joint <- predictions %>%
  mutate(specimenid = str_remove(sample_id,'^0+' # remove leading zeros
                                 ),
         specimenid = str_remove(specimenid,'..$' # remove last two characters
                                 )
         ) %>%
  # keep predictions, sample id, chorno age, and vars for filtering
  select(chrono_age, specimenid, age, imp_method, tissue) %>%
  # add in participant ids and participant info
  full_join(., 
            sample_long, 
            by = 'specimenid'
            )

# output 
################################################################################
write.csv(joint,
          file = 'data-processed/dnam-age-sample-info.csv',
          row.names = FALSE
          )