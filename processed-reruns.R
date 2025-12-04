## Author: Selene Banuelos
## Date: 11/25/2025
## Description: Identify timepoint 5 participants that have buccual telomere 
## data and those who were assayed more than once

# setup 
library(dplyr)
library(stringr)

# import data ##################################################################
# save raw data file names as list
file_names <- list.files(path = 'data-raw/t5/',
                         pattern = '.csv',
                         full.names = TRUE)

# import all raw data files and combine into one dataframe
raw_data <- purrr::map_dfr(file_names, 
                   # import .csv as df and add column with original file name
                   function(x) read.csv(x) %>% mutate(file_name = x)
                   )

# import master list of all T5 participants with buccal samples
all_t5 <- read.csv('data-raw/PEARLSBio-T5sWithBuccal_DATA_2025-12-02_1043.csv')

# clean up raw data ############################################################
clean_data <- raw_data %>%
  # create run_date column
  mutate(batch = str_extract(
    file_name, 
    '(?<=TEL_)(.*)(?= - Quant)') # (?<=prefix)(keep)(?=suffix)
    ) %>%
  select(c(Sample, Content, batch, contains('Cq'))) %>%
  # remove any rows corresponding to controls: Content = 'Std', 'Ctrl', 'NTC'
  filter(
    str_detect(Content, # column to filter on
                        "Std|Ctrl|NTC", 
                        negate = TRUE # keep rows that don't have these strings
                        )
  ) %>%
  # strip any leading zeros from sample names
  mutate(Sample = str_remove(
    Sample,
    '^0+') # match any number of zeros at the beginning of a string
    )

# identify reruns ##############################################################
reruns <- clean_data %>%
  group_by(Sample) %>%
  filter(n_distinct(batch) >1) %>% # keep samples that were run in > 1 batch
  distinct(Sample, batch)

# all samples run
run_ids <- unique(clean_data$Sample)

# all samples run more than once
rerun_ids <- unique(reruns$Sample)

# identify which participants were not assayed and why #########################
missing <- all_t5 %>%
  mutate(tel_data = case_when(specimenid %in% run_ids ~ 1,
                             .default = 0
                             )
         )

# output #######################################################################
write.csv(missing,
          'data-processed/missing-T5-tel.csv')