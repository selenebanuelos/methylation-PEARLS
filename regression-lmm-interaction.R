### Author: Selene Banuelos
### Date: 3/12/2026
### Description: Compare DNAm age deviation from chronological age across two
### timepoints between participants who experienced no PEARLS and those who
### experienced high PEARLS. Linear mixed effects models used to analyze
### this relationship in buccal and blood samples, separately. 

# setup
library(dplyr)
library(lme4)
library(purrr)

# import data 
################################################################################
# DNAm age data and other sample info
data <- read.csv('data-processed/dnam-age-sample-info.csv')

# demographics data (baseline (T2) household income: income_FPL_100)
demo <- read.csv('data-raw/pearls_dataset_2022-07-08.csv')

# data wrangling 
################################################################################
# combine all data needed
clean <- demo %>%
  filter(visitnum == 2) %>% # used only baseline data
  select(pearls_id, income_FPL_100) %>% # get household income data
  right_join(., data, by = 'pearls_id') %>% # combine all data together
  # make PEARLS categories no/high
  mutate(pearls = case_when(aces_baseline == 0 ~ 'no',
                            aces_baseline >= 5 ~ 'high'
                            )
         )

# factor all categorical variables
# timepoint
clean$timepoint <- factor(clean$timepoint,
                          levels = c(2, 5)) # reference: 2
# PEARLS status
clean$pearls <- factor(clean$pearls,
                       levels = c('no', 'high')) # reference: no
# sex
clean$sex <- factor(clean$sex, 
                    levels = c(0, 1),
                    labels = c('female', 'male')) # reference: female

# household income below 100% federal poverty level for household of 4 (<25k)
clean$income_FPL_100 <- factor(clean$income_FPL_100,
                               levels = c(0,1),
                               labels = c('no', 'yes')) # reference: no

# linear mixed effects models 
################################################################################
# blood samples
blood_pedbe <- lmer(ped_be ~ age*pearls + pearls + age + sex  + income_FPL_100 + (1 | pearls_id),
                     data = filter(clean, tissue == 'blood'))

blood_horvath2 <- lmer(horvath2 ~ age*pearls + pearls + age + sex + income_FPL_100 + (1 | pearls_id),
                     data = filter(clean, tissue == 'blood'))

blood_horvath1 <- lmer(horvath_age ~ age*pearls + pearls + age + sex + income_FPL_100 + (1 | pearls_id),
                          data = filter(clean, tissue == 'blood'))

# buccal samples
buccal_pedbe <- lmer(ped_be ~ age*pearls + pearls + age + sex + income_FPL_100 + (1 | pearls_id),
                     data = filter(clean, tissue == 'buccal'))

buccal_horvath2 <- lmer(horvath2 ~ age*pearls + pearls + age + sex + income_FPL_100 + (1 | pearls_id),
                       data = filter(clean, tissue == 'buccal'))

buccal_horvath1 <- lmer(horvath_age ~ age*pearls + pearls + age + sex + income_FPL_100 + (1 | pearls_id),
                       data = filter(clean, tissue == 'buccal'))

# create table of estimates (95% CI) for all terms in all models
################################################################################
# because of small sample size, apply Kenward-Roger approximation to the degrees 
# of freedom to calculate 95% CI
get_results <- function(model, # model object
                        tissue # string specifying tissue type
                        ) {
  
# Kenward-Roger approximation for 95% CI
cis <- parameters::ci_kenward(model) %>%
  # round CI bounds to 2 decimal places
  mutate(across(c('CI_low', 'CI_high'), 
                \(x) round(x, digits = 2)
                )
         )

# get coefficient estimates from model object
estimates <- fixef(model) %>% # produces named numeric vector
  # convert to data frame to facilitate joining
  data.frame(.) %>% 
  # convert row names into columns to facilitate joining
  tibble::rownames_to_column(., 'Parameter') %>%
  # rename column that contains coefficient estimates
  rename(beta = '.') %>%
  # round coefficient estimates to 2 decimal places
  mutate(beta = round(beta, digits = 2))

# get response variable name
y <- names(model@frame)[1]

# join coefficient estimates and corresponding 95% CIs into one table
full_join(estimates, cis, by = 'Parameter') %>%
  # create column containing estimate + full CI from bounds
  mutate(beta_95_ci = paste0(beta, 
                             ' (', 
                             CI_low, 
                             ', ', 
                             CI_high, 
                             ')'
                             )
         ) %>%
  # keep only desired variables for table
  select(Parameter, beta_95_ci) %>%
  # add in response variable name
  mutate(outcome = y,
         tissue = tissue)
}

# put all model objects into list
blood <- list(blood_horvath1, blood_horvath2, blood_pedbe) %>%
  # create data frame with combined results from all models
  map_df(., get_results, 'blood')

buccal <- list(buccal_horvath1, buccal_horvath2, buccal_pedbe) %>%
  map_df(., get_results, 'buccal')

# put all results into one data frame
results <- rbind(blood, buccal)

# output
################################################################################
write.csv(results, 'data-processed/mm-interaction-results.csv', row.names = FALSE)