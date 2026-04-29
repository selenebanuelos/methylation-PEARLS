### Author: Selene Banuelos
### Date: 4/29/2026
### Description: Make table of PEARLS items reported by pilot study participants

# setup
library(dplyr)

# import data
################################################################################
# dataset that includes answers to PEARLS screener
all_data <- read.csv('data-raw/pearls_dataset_2022-07-08.csv')

# sample info from participants in pilot study
participants <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# data wrangling
################################################################################
# make vector of original 10 ACEs items
aces <- c('incarceration',
          'emotional_abuse',
          'mental_illness',
          'verbal_abuse',
          'substance_abuse',
          'physical_neglect',
          'divorce_cohesion',
          'domestic_violence',
          'physical_abuse',
          'sexual_abuse')

# make vector of additional 7 related life events
rle <- c('neighborhood_violence', 
         'food_insecurity',
         'discrimination',
         'housing_insecurity',
         'physical_illness',
         'separation_foster_immi',
         'caregiver_death')

# make vector of pilot study participant IDs with reported ACEs
reported_aces <- filter(participants, aces_baseline != 0) %>%
  .$subjectid

# all PEARLS screener data from baseline
pearls <- all_data %>%
  # keep ID, visit number, PEARLS items
  select(pearls_id, visitnum, contains(aces), contains(rle)) %>%
  # only look at participants from pilot study that reported ACEs
  filter(pearls_id %in% reported_aces,
         # interested in baseline ACEs only
         visitnum == 1 | visitnum == 2)

# check if answers were collected at initial recruitment or biospecimen visit
visit_check <- pearls %>%
  # count number of PEARLS items reported at each visit (there are repeated variables)
  mutate(total_pearls = rowSums(.[,-c(1,2)], na.rm = TRUE)) %>%
  relocate(total_pearls, .after = visitnum)
# looks like all PEARLS screeners were conducted during visit number 1 (recruitment)

# check for differences between different versions of PEARLS variables
var_versions <- pearls %>%
  # PEARLS were evaluated at first visit
  filter(visitnum == 1)
# it looks like the "ident_plus_deident" versions contain complete data

# cleaned up dataset with all PEARLS items
clean_pearls <- pearls %>%
  filter(visitnum == 1) %>%
  select(pearls_id, visitnum, contains('ident_plus_deident')) %>%
  mutate_if(is.numeric, as.factor)

# remove "_ident_plus_deident" suffix from variable names
names(clean_pearls) <- sub('_ident_plus_deident', '', names(clean_pearls))

# data visualization
################################################################################
# make dotplot or heatmap to show how many participants reported which events
