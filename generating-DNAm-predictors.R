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

# data wrangling ###############################################################
# clean up phenotype data for use in methscore() below
blood_pheno <- blood_info %>%
  # originally, 0=Female. Change to 1=Female, 0=Male
  mutate(Female = case_when(sex_.0.F. == '0' ~ 1,
                            sex_.0.F. == '1' ~ 0)
         ) %>%
  dplyr::rename('SampleID' = 'Sample_Name',
         'Age' = 'age_baseline'
         ) %>%
  select(c(SampleID, Age, Female))

buccal_pheno <- buccal_info %>%
  # originally, 0=Female. Change to 1=Female, 0=Male
  mutate(Female = case_when(sex_.0.F. == '0' ~ 1,
                            sex_.0.F. == '1' ~ 0)
  ) %>%
  dplyr::rename('SampleID' = 'Sample_Name',
         'Age' = 'age_baseline'
  ) %>%
  select(c(SampleID, Age, Female))
 
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