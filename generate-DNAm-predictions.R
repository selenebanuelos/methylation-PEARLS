### Author: Selene Banuelos
### Date: 11/25/2025
### Description: Generate methylation-based predictions for participants

# setup
library(dplyr)
library(ENmix) #BiocManager::install('ENmix')

# import data ##################################################################
# normalized beta values
blood <- readRDS('data-raw/Final_Cleaned_Blood_Betas_n39.rds')
buccal <- readRDS('data-raw/Final_Cleaned_Buccal_Betas_n38.rds')

# import processed sample info
blood_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')
buccal_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# import participants demographics data (contains chronological age)
demo <- read.csv('data-raw/pearls_data_LauraDiaz_2025_11_20.csv')

# data wrangling ###############################################################
# clean up chronological age data
age <- demo %>%
  select(pearls_id, contains('collectionage')) %>%
  tidyr::pivot_longer(
    cols = !pearls_id,
    names_to = c('visitnum'),
    names_pattern = '_t(.)', # keep number after '_t', which specifies timepoint
    values_to = 'Age'
  ) %>%
  # get rid of t4 info, we don't have data for that timepoint
  filter(visitnum != 4) 
  

# clean up phenotype data for use in methscore() below
# put blood and buccal sample info into list to clean up
sample_info <- list(blood_info, buccal_info)

# vector of variables needed to generate predictions
vars <- c('subjectid', 'Timepoint', 'sex_.0.F.', 'Sample_Name', 'Tissue')

# clean up phenotype data
phenotype <- sample_info %>%
  # select only relevant variables and combine into one df
  purrr::map_df(.,
                function(x) x %>% select(vars)
  ) %>%
  # originally, 0=Female. Change to 1=Female, 0=Male
  mutate(Female = case_when(sex_.0.F. == '0' ~ 1,
                            sex_.0.F. == '1' ~ 0)
  ) %>%
  # rename for joining
  dplyr::rename('pearls_id' = 'subjectid') %>%
  # remove 'T' from timepoint value for joining
  mutate(visitnum = str_remove(Timepoint, '^T')
  ) %>%
  # add in age at biospecimen collection
  left_join(.,
            age,
            by = c('pearls_id', 'visitnum')
            ) %>%
  # rename for methscore() below
  dplyr::rename('SampleID' = 'Sample_Name')

blood_pheno <- phenotype %>%
  filter(Tissue == 'Blood') %>%
  # keep only relevant variables for methscore() below
  select(SampleID, Age, Female)

buccal_pheno <- phenotype %>%
  filter(Tissue == 'Buccal') %>%
  # keep only relevant variables for methscore() below
  select(SampleID, Age, Female)
 
# predict ######################################################################
# predict DNA methylation age & plasma protein levels

# using mean methylation values to impute missing CpGs
blood_pred_mean <- methscore(datMeth = blood, # methylation beta value matrix
                              datPheno = blood_pheno, # phenotype data,
                              fastImputation = TRUE, # use mean meth values
                              normalize = FALSE # data previously normalized
                              ) %>%
  mutate(tissue = 'blood',
         imp_method = 'mean')

buccal_pred_mean <- methscore(datMeth = buccal,
                              datPheno = buccal_pheno,
                              fastImputation = TRUE, # use mean meth values
                              normalize = FALSE 
                              ) %>%
  mutate(tissue = 'buccal',
         imp_method = 'mean')

# using KNN to impute missing CpGs
blood_pred_knn <- methscore(datMeth = blood, 
                             datPheno = blood_pheno, 
                             fastImputation = FALSE, # use KNN imputation
                             normalize = FALSE 
                            ) %>%
  mutate(tissue = 'blood',
         imp_method = 'knn')

buccal_pred_knn <- methscore(datMeth = buccal,
                              datPheno = buccal_pheno,
                              fastImputation = FALSE, # use KNN imputation
                              normalize = FALSE 
                             ) %>%
  mutate(tissue = 'buccal',
         imp_method = 'knn')

# output #######################################################################
# combine into one dataframe
dnam_predictions <- rbind(blood_pred_mean, 
                         buccal_pred_mean, 
                         blood_pred_knn, 
                         buccal_pred_knn)

# save as csv
write.csv(dnam_predictions,
          'data-processed/methylation-predictions.csv',
          row.names = FALSE)