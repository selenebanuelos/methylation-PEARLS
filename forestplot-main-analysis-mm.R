### Author: Selene Banuelos
### Date: 5/27/2026
### Description: Forest plot to visualize regression results using mixed models

# setup
library(dplyr)
library(tidyr)
library(ggplot2)

# import data
################################################################################
# import regression results
results <- read.csv('data-processed/mm_results.csv')

# data wrangling
################################################################################
# split estimate and 95% CI bounds into separate columns
split <- results %>%
  extract(beta_95_ci, # split this column 
          # into these 3 columns:
          into = c('estimate', 'ci_low', 'ci_high'),
          # keep groups: estimate, lower bound, upper bound
          regex = '([\\d.-]+)\\s*\\(([\\d.-]+),\\s*([\\d.-]+)\\)'
           ) %>%
  mutate_at(c('estimate', 'ci_low', 'ci_high'), as.numeric)

# control order and labels of factored variables
# covariates
split$Parameter <- factor(split$Parameter, 
                          levels = c('(Intercept)',
                                     'income_FPL_100yes',
                                     'sexmale',
                                     'age:pearlshigh',
                                     'pearlshigh',
                                     'age'),
                          labels = c('Intercept',
                                     'Household Income',
                                     'Sex',
                                     'Chrono Age*PEARLS',
                                     'PEARLS',
                                     'Chrono Age'))
# clocks
split$outcome <- factor(split$outcome, 
                        levels = c('ped_be',
                                   'horvath2',
                                   'horvath_age'),
                        labels = c('PedBE',
                                   'Skin & Blood',
                                   'Multi-tissue'))
# tissue type
split$tissue <- factor(split$tissue,
                       levels = c('buccal',
                                  'blood'),
                       labels = c('Buccal',
                                  'Blood'))
# data visualization
################################################################################
# forest plot
plot <- split %>%
  # don't include intercept in plot
  filter(Parameter != 'Intercept') %>%
  ggplot(aes(x = Parameter, # show point est & 95% CI for each coefficient
             y = estimate,
             ymin = ci_low,
             ymax = ci_high)
         ) +
  # plot points with bars to represent CI
  geom_pointrange(color = 'black', size = 0.5) +
  # plot line at null value of 0
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  # show coefficient names on y-axis and values on x-axis
  coord_flip() +
  # no axis label for coefficient names
  xlab(" ") +
  ylab("Coefficient Estimate (95% Confidence Interval)") +
  # separate plots for tissue/clock combinations
  facet_grid(tissue ~ outcome, switch = 'y') + 
  # formatting
  theme_light() +
  theme(panel.border = element_rect(color = "#CCCCCC", fill = NA, linewidth = 2),
        strip.text = element_text(face = 'bold', size = 14),
        strip.background = element_rect(fill = 'darkgrey'),
        axis.text = element_text(size = 12)
        )

# output
################################################################################
ggsave('figures/forestplot-mm.png', 
       plot,
       width = 9,
       height = 5.5,
       units = 'in')