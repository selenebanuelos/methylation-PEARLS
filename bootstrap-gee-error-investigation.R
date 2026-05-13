### Author: Selene Banuelos
### Date: 4/30/2026
### Description: investigate errors and warnings produced by model fits

# setup
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
options(scipen = 999)

# import data 
################################################################################
# load bootstrap estimates and observed model fits
load('data-processed/bootstrap-estimates.RData')

# check boostrap distributions for skew
################################################################################
check_bias <- function(boot_df, # df of bootstrapped estimates
                       obs_fit # model object fit with observed data
){
  # boot_df = blood_boot
  # obs_fit = obs_blood_sb
  
  
  # calculate bias: mean(bootstrap estimates) - observed estimate
  point <- round(
    obs_fit$coefficients["pearlshigh"], 
    digits = 3
  )
  cat("Point estimate:", point, "\n")
  
  boot_mean <- round(
    mean(boot_df$pearlshigh, na.rm = TRUE),
    digits = 3
  )
  cat("Bias (mean bootsrap estimate - observed point estimate):", boot_mean - point, "\n")
  cat("Skewness: sign tells you which tail is longer\n")
  
  # visualize bootstrap coefficient estimate distributions
  boot_df %>%
    # make data longer for plotting
    pivot_longer(
      cols = -c(error, warning),
      names_to = "Covariate",
      values_to = "Estimate"
    ) %>%
    # plot distributions of all covariates
    ggplot(aes(x = Estimate)) +
    geom_density() +
    facet_wrap(~Covariate)
}

# blood samples
check_bias(blood_sb_boot, obs_blood_sb)
check_bias(blood_pbe_boot, obs_blood_pbe)

# buccal samples
check_bias(buccal_sb_boot, obs_buccal_sb)
check_bias(buccal_pbe_boot, obs_buccal_pbe)
# little bit of bias in both bootstrapped distributions of estimate

# check proportion of bootstrap models generated warning messages and errors
################################################################################
# count number of bootstrap models with and without warnings about model fit
# in blood samples
count_warning <- function(boot_df){
  
  # get name of dataframe, which specifies sample type and epigenetic clock used
  data_name <- as.character(substitute(boot_df))
  
  # create new variable name for proportion data that included sample/clock info
  sample_prop <- paste0(data_name, '_warning_prop')
  
  # return dataframe with error type counts
  boot_df %>% 
    count(warning) %>% 
    # calculate proportions of errors from all resamples
    mutate(!!sample_prop := n/10000) %>%
    # rename 
    rename(!!data_name := n)
}

count_error <- function(boot_df){
  
  # get name of dataframe, which specifies sample type and epigenetic clock used
  data_name <- as.character(substitute(boot_df))
  
  # create new variable name for proportion data that included sample/clock info
  sample_prop <- paste0(data_name, '_error_prop')
  
  # return dataframe with error type counts
  boot_df %>% 
    count(error) %>% 
    # calculate proportions of errors from all resamples
    mutate(!!sample_prop := n/10000) %>%
    # rename 
    rename(!!data_name := n)
}

# get proportions of different types of warnings with model fits in bootstrapping
warnings <- list(count_warning(blood_sb_boot),
                 count_warning(blood_pbe_boot),
                 count_warning(buccal_sb_boot),
                 count_warning(buccal_pbe_boot)
                 ) %>%
  reduce(full_join, by = "warning")

# get proportions of different types of errors with model fits in bootstrapping
errors <- list(count_error(blood_sb_boot),
               count_error(blood_pbe_boot),
               count_error(buccal_sb_boot),
               count_error(buccal_pbe_boot)
               ) %>%
  reduce(full_join, by = "error")

# calculate 95% CI after removing estimates from models that generated warnings
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

# build quantile 95% CI from bootstrap estimates from model fits with no warnings
blood_sb_boot %>%
  filter(warning == 'none') %>%
  quantile_ci()
  

blood_pbe_boot %>%
  filter(warning == 'none') %>%
  quantile_ci()

buccal_sb_boot %>%
  filter(warning == 'none') %>%
  quantile_ci()

buccal_pbe_boot %>%
  filter(warning == 'none') %>%
  quantile_ci()

# output 
################################################################################
# save warning counts and proportions
write.csv(warnings, 
          file = 'data-processed/bootstrap-warnings.csv', 
          row.names = FALSE)

# save error counts and proportions
write.csv(errors, 
          file = 'data-processed/bootstrap-errors.csv', 
          row.names = FALSE)