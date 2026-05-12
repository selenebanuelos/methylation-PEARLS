### Author: Selene Banuelos
### Date: 4/6/2026
### Description: calculate mean, min, and max of time between PEALRS visits

library(dplyr)
library(tidyr)

data <- read.csv('data-processed/analysis-ready-dataset.csv')

# difference in age between two timepoints blood samples
diff_blood <- data %>% 
  filter(tissue == 'blood') %>%
  pivot_wider(id_cols = pearls_id,
              names_from = timepoint,
              names_glue = 'T{timepoint}',
              values_from = age) %>%
  mutate(diff = T5 - T2)

mean(diff_blood$diff, na.rm = TRUE)
min(diff_blood$diff, na.rm = TRUE)
max(diff_blood$diff, na.rm = TRUE)

# difference in age between two timepoints blood samples
diff_buccal <- data %>% 
  filter(tissue == 'buccal') %>%
  pivot_wider(id_cols = pearls_id,
              names_from = timepoint,
              names_glue = 'T{timepoint}',
              values_from = age) %>%
  mutate(diff = T5 - T2)

mean(diff_buccal$diff, na.rm = TRUE)
min(diff_buccal$diff, na.rm = TRUE)
max(diff_buccal$diff, na.rm = TRUE)