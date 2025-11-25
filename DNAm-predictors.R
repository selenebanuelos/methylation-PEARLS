### Author: Selene Banuelos
### Date: 11/25/2025
### Description: Generate methylation predictors for participants

# setup
library(dplyr)
library(ENmix)

# import data
blood <- readRDS('data-raw/Final_Cleaned_Blood_Betas_n39.rds')
buccal <- readRDS('data-raw/Final_Cleaned_Buccal_Betas_n38.rds')

# import processed sample info
blood_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')
buccal_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# clean up phenotype data for use in methscore() below
blood_pheno <- blood_info %>%
  # originally, 0=Female. Change to 1=Female, 0=Male
  mutate(Female = case_when(sex_.0.F. == '0' ~ 1,
                            sex_.0.F. == '1' ~ 0)
         ) %>%
  rename('SampleID' = 'Sample_Name',
         'Age' = 'age_baseline'
         ) %>%
  select(c(SampleID, Age, Female))

buccal_pheno <- buccal_info %>%
  # originally, 0=Female. Change to 1=Female, 0=Male
  mutate(Female = case_when(sex_.0.F. == '0' ~ 1,
                            sex_.0.F. == '1' ~ 0)
  ) %>%
  rename('SampleID' = 'Sample_Name',
         'Age' = 'age_baseline'
  ) %>%
  select(c(SampleID, Age, Female))
 
# calculate methylation predictors: DNA methylation age & plasma protein levels
blood_predictors <- methscore(datMeth = blood, # methylation beta value matrix
                              datPheno = blood_pheno, # phenotype data
                              normalize = FALSE # data previously normalized
                              )

buccal_predictors <- methscore(datMeth = buccal,
                              datPheno = buccal_pheno,
                              normalize = FALSE # data previously normalized
)