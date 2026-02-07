### Author: Selene Banuelos
### Date: 11/25/2025
### Description: Generate methylation-based predictions for participants for 
### each timepoint, independently

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
  mutate(visitnum = stringr::str_remove(Timepoint, '^T')
  ) %>%
  # add in age at biospecimen collection
  left_join(.,
            age,
            by = c('pearls_id', 'visitnum')
            ) %>%
  # rename for methscore() below
  dplyr::rename('SampleID' = 'Sample_Name')

# separate phenotype data by timepoint and tissue type
blood_pheno_2 <- phenotype %>%
  filter(Tissue == 'Blood',
         Timepoint == 'T2'
         ) %>%
  # keep only relevant variables for methscore() below
  select(SampleID, Age, Female)

blood_pheno_5 <- phenotype %>%
  filter(Tissue == 'Blood',
         Timepoint == 'T5'
         ) %>%
  select(SampleID, Age, Female)

buccal_pheno_2 <- phenotype %>%
  filter(Tissue == 'Buccal',
         Timepoint == 'T2'
         ) %>%
  select(SampleID, Age, Female)

buccal_pheno_5 <- phenotype %>%
  filter(Tissue == 'Buccal',
         Timepoint == 'T5'
         ) %>%
  select(SampleID, Age, Female)

# separate methylation data by timepoint and tissue type
blood_2 <- blood[, blood_pheno_2$SampleID]
blood_5 <- blood[, blood_pheno_5$SampleID]
buccal_2 <- buccal[, buccal_pheno_2$SampleID]
buccal_5 <- buccal[, buccal_pheno_5$SampleID]
 
# predict ######################################################################
# predict epigenetic age & plasma protein levels
# using mean methylation values to impute missing CpGs

# blood samples from t2
blood_pred_2 <- methscore(datMeth = blood_2, # methylation beta value matrix
                          datPheno = blood_pheno_2, # phenotype data,
                          fastImputation = TRUE, # use mean meth values
                          normalize = FALSE # data previously normalized
                          ) %>%
  mutate(tissue = 'blood', 
         imp_method = 'mean'
         )

# blood samples from t5
blood_pred_5 <- methscore(datMeth = blood_5, 
                          datPheno = blood_pheno_5, 
                          fastImputation = TRUE, 
                          normalize = FALSE 
                          ) %>%
  mutate(tissue = 'blood', 
         imp_method = 'mean'
         )

# buccal samples from t2
buccal_pred_2 <- methscore(datMeth = buccal_2,
                           datPheno = buccal_pheno_2,
                           fastImputation = TRUE, 
                           normalize = FALSE 
                           ) %>%
  mutate(tissue = 'buccal',
         imp_method = 'mean'
         )

# buccal samples from t5
buccal_pred_5 <- methscore(datMeth = buccal_5,
                           datPheno = buccal_pheno_5,
                           fastImputation = TRUE, 
                           normalize = FALSE 
                           ) %>%
  mutate(tissue = 'buccal',
         imp_method = 'mean'
         )

# output #######################################################################
# combine into one dataframe
dnam_predictions <- rbind(blood_pred_2,
                          blood_pred_5,
                          buccal_pred_2,
                          buccal_pred_5
                          )

# save as csv
write.csv(dnam_predictions,
          'data-processed/methylation-predictions.csv',
          row.names = FALSE)