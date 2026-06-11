### Author: Selene Banuelos
### Date: 6/10/2026
### Description: Investigate PEARLS domains in the 10 high adversity participants

# setup
library(dplyr)
library(tidyr)
library(ggplot2)

# import data
################################################################################
# PEARLS 3 domains data
domains <- read.csv('data-raw/pearls_domains_SeleneB_2026-06-01.csv')

# pilot study participants
participants <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# data wrangling
################################################################################
# get baseline PEARLS domains for pilot study participants
baseline <- domains %>%
  # used PEARLS at baseline as exposure in pilot study
  filter(visitnum == 1) %>%
  # rename for joining
  rename(subjectid = pearls_id) %>%
  # remove unwanted variables
  select(-visitnum) %>%
  # join domain data to pilot study participant data
  right_join(., participants, by = 'subjectid') %>%
  # keep only participants who reported PEARLS items
  filter(aces_baseline >= 5) %>%
  # make indicator variable for each PEARLS domain
  mutate(maltreatment_bin = ifelse(maltreatment > 0, 1, 0),
         household_chal_bin = ifelse(household_chal > 0, 1, 0),
         sdoh_bin = ifelse(sdoh > 0, 1, 0)
         )
# 9/10 of participants reported at least 1 item in all 3 domains

# make data longer for plotting
long <- baseline %>%
  select(subjectid, maltreatment_bin, household_chal_bin, sdoh_bin) %>%
  pivot_longer(cols = -subjectid,
               names_to = 'domain',
               values_to = 'presence'
  ) %>%
  # factor indicators of reporting each PEARLS domain
  mutate(presence = factor(presence, 
                           levels = c(0,1),
                           labels = c('Absent',
                                      'Present'))) %>%
  group_by(subjectid) %>%
  # mask participant ID with index value
  mutate(participant_index = factor(cur_group_id())) %>%
  ungroup()

# data visualization
################################################################################
# make dotplot to show which participants reported which events
plot <- ggplot(long, aes(
  x = participant_index,
  y = domain,
  # fill in dots depending on if domain was reported or not
  fill = presence)
  ) +
  geom_point(shape = 21, # circle
             size = 7, # size of circles
             color = '#BCE9C5FF', # outline of circles
             stroke = 1.5 # width of outline
  ) +
  # change colors of dots/outlines
  scale_fill_manual(values = c('Absent' = 'white',
                               'Present' = '#18BDB0FF')
  ) +
  # formatting
  theme_minimal() +
  labs(
    title = 'PEARLS domains reported by each participant',
    x = 'Participant',
    y = 'Domain',
    fill = 'Presence' 
  ) +
  theme(panel.grid = element_blank()) +
  # format domain labels
  scale_y_discrete(labels = c('sdoh_bin' = 'Social Context', 
                              'maltreatment_bin' = 'Maltreatment',
                              'household_chal_bin' = 'Household Challenges')
                   )

# output
################################################################################
# save plot as png
ggsave('figures/dotplot-pearls-domains.png', 
       plot,
       width = 8,
       height = 2,
       units = 'in')