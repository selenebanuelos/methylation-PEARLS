### Author: Selene Banuelos
### Date: 4/29/2026
### Description: Make heatmap or dotplot type plot to show prevalence of each PEARLS item
### in sample participants that reported >= 5 PEARLS
### Notes: https://medium.com/@one-more-step/how-to-draw-an-ordered-binary-presence-absence-dot-plot-or-heat-map-in-r-512feebd19d8

# setup
library(dplyr)
library(tidyr)
library(ggplot2)

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
  select(pearls_id, contains('ident_plus_deident')) %>%
  mutate_if(is.numeric, as.factor)

# remove "_ident_plus_deident" suffix from variable names
names(clean_pearls) <- sub('_ident_plus_deident', '', names(clean_pearls))

# make data longer for plotting
long <- clean_pearls %>%
  pivot_longer(cols = -pearls_id,
               names_to = 'item',
               values_to = 'presence'
  ) %>%
  mutate(presence = factor(presence, 
                           levels = c(0,1),
                           labels = c('Absent',
                                      'Present'))) %>%
  group_by(pearls_id) %>%
  mutate(pearls_score = sum(presence == 'Present'),
         participant_index = factor(cur_group_id())) %>%
  ungroup() %>%
  group_by(item) %>%
  mutate(item_freq = sum(presence == 'Present')) %>%
  ungroup()

# data visualization
################################################################################
# make dotplot to show which participants reported which events
plot <- ggplot(long, aes(
  # order participants by PEARLS score, highest on left and descending to right
  x = reorder(participant_index, pearls_score, decreasing = TRUE),
  # order items with highest reported freq on top and descending towards bottom
  y = reorder(item, item_freq),
  # fill in dots depending on if item was reported or not
  fill = presence)
  ) +
  geom_point(shape = 21, # circle
             size = 7, # size of circles
             color = '#BCE9C5FF', # outline of circles
             stroke = 1.5 # width of outline
             ) +
  scale_fill_manual(values = c('Absent' = 'white',
                               'Present' = '#18BDB0FF')
                    ) +
  theme_minimal() +
  labs(
    title = 'PEARLS items reported by each participant',
    x = 'Participant',
    y = 'Item',
    fill = 'Presence' 
    ) +
  theme(panel.grid = element_blank())

# output
################################################################################
ggsave('figures/dotplot-pearls-items.png', 
       plot,
       width = 6.5,
       height = 6,
       units = 'in')