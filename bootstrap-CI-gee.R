### Author: Selene Banuelos
### Date: 4/5/2026
### Description: Generate quantile-based 95% CI from bootstrap estimates
### some useful info: https://r-statistics.co/Bootstrap-Confidence-Intervals-in-R.html

# setup
library(dplyr)
options(scipen = 999)

# import data 
################################################################################
# load data frames with bootstrap estimates for blood and buccal samples
load('data-processed/bootstrap-estimates.RData')

# data wrangling
################################################################################
# multiply bootstrap estimates by 5 to generate CI for parameter of interest
blood_sb_param <- blood_sb_boot %>%
  select('age:pearlshigh', error, warning) %>%
  mutate(times_5 = `age:pearlshigh`*5)

# calculate 95% confidence intervals using 2.5th and 97.5th quantiles
################################################################################
# function that builds qunatile-based 95% CI from bootstrap distributions
quantile_ci <- function(boot_df) {
  
  # format bootstrap estimates
  boot_est <- boot_df %>%
    # remove rows for models that could not be fit (errored out)
    filter(error == 'none') %>%
    # remove error and warning columns
    select(-c(error, warning))
  
  # return data frame containing  95% CI bounds for each estimate
  apply(X = boot_est, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)) %>%
    round(., digits = 2) 
}

# build quantile 95% confidence intervals for estimates of all coefficients
blood_sb_boot %>%
  # remove estimates from model fits that generated warnings?
  quantile_ci()

blood_pbe_boot %>%
  quantile_ci()

buccal_sb_boot %>%
  quantile_ci()

buccal_pbe_boot %>%
  quantile_ci()

# build quantile 95% confidence intervals for parameter of interest
blood_sb_param %>%
  quantile_ci()

# calculate Wald-type 95% confidence intervals 
################################################################################
# function that creates 95% CI using robust SE
wald_ci <- function(boot_df){
  
  # format bootstrap estimates
  boot_est <- boot_df %>%
    # remove rows for models that could not be fit (errored out)
    filter(error == 'none') %>%
    # remove error and warning columns
    select(-c(error, warning))
  
  # sample size
  n <- nrow(boot_df)
  
  # calculate bounds
  lower_bound <- function(x) {mean(x) - 1.96 * sd(x)}
  upper_bound <- function(x) {mean(x) + 1.96 * sd(x)}
  
  # create dataframe with 95% CI bounds
  boot_est %>%
    summarise(across(everything(), list(lower = lower_bound, 
                                        upper = upper_bound
                                        )
                     )
              ) %>%
    round(., digits = 2)
}

# build Wald-type 95% CI for all coefficients
wald_ci(blood_sb_boot)
wald_ci(blood_pbe_boot)
wald_ci(buccal_sb_boot)
wald_ci(buccal_pbe_boot)

# build Wald-type 95% CI for parameter of interest
wald_ci(blood_sb_param)