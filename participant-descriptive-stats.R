### Author: Selene Banuelos
### Date: 1/14/2025
### Description: Summarize participant characteristics by adversity status

# setup
library(dplyr)
library(table1)

# import data
################################################################################
# sample information
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# data on partipant chronological age at each timepoint
age <- read.csv('data-raw/pearls_data_LauraDiaz_2025_11_20.csv')

# data on caregiver education at baseline
caregiver_edu <- read.csv('data-raw/pearls_dataset_2022-07-08.csv') %>%
  filter(visitnum == 1) %>%
  select(pearls_id, caregiver_edu_4groups)

# samples that passed qc
passed_qc <- read.csv('data-processed/samples-passed-qc.csv') %>%
  select(specimenid, passed_qc, tissue) %>%
  mutate(specimenid = as.character(specimenid))
  
# data wrangling
################################################################################
# reformat age data longer to separate age at biospecimen collection & timepoint
age_long <- age %>%
  select(pearls_id, contains('collectionage_')) %>%
  tidyr::pivot_longer(
    cols = !pearls_id,
    names_pattern = 'collectionage_t(.)', # only keep timepoint number 
    names_to = 'timepoint',
    values_to = 'age'
  )

# clean up DNAme sample info with specimen id, sex, ACEs at baseline
sample_long <- sample %>%
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
               names_pattern = 'T(.)_(.*)', # don't keep 'T' prefix for timepoint
               names_to = c('timepoint', '.value')
  ) %>%
  mutate(specimenid = as.character(specimenid)) %>%
  dplyr::rename(pearls_id = subjectid) %>%
  # create categorical ACEs var (no ACEs/high ACEs) from ACE score at BASELINE
  mutate(aces_cat = case_when(aces_baseline == 0 ~ "none",
                              aces_baseline != 0 ~ "high"
                              )
         )

# combine all variables into one dataframe
# start with specimend id, sex, ACEs at baseline
characteristics <- select(sample_long, !age_baseline) %>%
  # add in chronological age at each timepoint
  left_join(age_long,
            by = c('pearls_id', 'timepoint')
            ) %>%
  # add in caregiver education at baseline
  left_join(caregiver_edu,
            by = 'pearls_id'
            ) %>%
  # add in indicator of sample passing QC
  left_join(passed_qc,
            by = 'specimenid'
            )
# will need to add in cell type eventually? but that might take some extra work
# will need to review email from Kobor lab b/c I think they mentioned that cell
# type proportions should be estimated after batch correction? not sure...
# cell type proportions are reported in DNAme sample informationl procided by 
# Kobor lab

# data visualization
################################################################################
# create descriptive statistics tables (table 1), stratified by adversity status
# table for blood samples
stats_blood <- table1(
  # display the following characteristics, stratified on ACEs status
  ~ factor(sex) + age + caregiver_edu_4groups | aces_cat, 
  # specify which participants to include in table
  data = characteristics %>%
    # keep only blood data
    filter(tissue == 'blood') %>%
    # only include participants that passed QC
    filter(passed_qc == 1),
  caption = 'Descriptive statistics for blood samples'
  )

# table for buccal samples
stats_buccal <- table1(
  # display the following characteristics, stratified on ACEs status
  ~ factor(sex) + age + caregiver_edu_4groups | aces_cat, 
  # specify which participants to include in table
  data = characteristics %>%
    # keep only blood data
    filter(tissue == 'buccal') %>%
    # only include participants that passed QC
    filter(passed_qc == 1),
  caption = 'Descriptive statistics for buccal samples'
)

# output
