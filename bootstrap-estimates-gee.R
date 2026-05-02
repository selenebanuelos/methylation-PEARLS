### Author: Selene Banuelos
### Date: 4/5/2026
### Description: Generate bootstrap estimates of regression coefficients
### from marginal model (GEE)

# setup
library(dplyr)
library(tidyr)
library(gee)
library(ggplot2)
options(scipen = 999)

# import data 
################################################################################
# analysis ready dataset
data <- read.csv('data-processed/analysis-ready-dataset.csv')

# data wrangling
################################################################################
analysis_vars <- c('pearls_id', 
                   'income_FPL_100', 
                   'horvath2', 
                   'ped_be', 
                   'age', 
                   'tissue', 
                   'sex', 
                   'timepoint', 
                   'pearls')

# factor variables
# pearls ID: needs to be factored to use with gee()
data$pearls_id <- factor(data$pearls_id)

# timepoint
data$timepoint <- factor(data$timepoint,
                          levels = c(2, 5)) # reference: 2
# PEARLS status
data$pearls <- factor(data$pearls,
                       levels = c('no', 'high')) # reference: no
# sex
data$sex <- factor(data$sex, 
                    levels = c('female', 'male')) # reference: female

# household income below 100% federal poverty level for household of 4 (<25k)
data$income_FPL_100 <- factor(data$income_FPL_100,
                               levels = c('no', 'yes')) # reference: no

# make data wider so we can do cluster-based bootstrap sampling (sample both
# timepoints per individual)
wide_data <- select(data, all_of(analysis_vars)) %>%
  pivot_wider(
    id_cols = c(pearls_id, tissue), 
    names_from = timepoint,
    values_from = c(income_FPL_100, horvath2, ped_be, age, sex, pearls),
    names_glue = "{.value}_{timepoint}"
  )
    

# bootstrapping setup
################################################################################
# fit model with observed data
obs_blood_fit <- gee(horvath2 ~ age*pearls + pearls + age + sex + income_FPL_100,
              id = pearls_id,
              data = filter(data, tissue == 'blood'),
              family = 'gaussian',
              corstr = 'exchangeable')

obs_buccal_fit <- gee(horvath2 ~ age*pearls + pearls + age + sex + income_FPL_100,
                     id = pearls_id,
                     data = filter(data, tissue == 'buccal'),
                     family = 'gaussian',
                     corstr = 'exchangeable')

# split data into no/high adversity datasets to do stratified boostrap sampling
################################################################################
# sample stratified for no PEARLS participants
blood_no <- wide_data %>%
  # sample from participants with PEARLS = 0
  filter(tissue == 'blood',
         pearls_2 == 'no' | pearls_5 == 'no')

# sample stratified for high PEARLS participants
blood_high <- wide_data %>%
  filter(tissue == 'blood',
         pearls_2 == 'high' | pearls_5 == 'high')

# sample stratified for no PEARLS participants
buccal_no <- wide_data %>%
  # sample from participants with PEARLS = 0
  filter(tissue == 'buccal',
         pearls_2 == 'no' | pearls_5 == 'no')

# sample stratified for high PEARLS participants
buccal_high <- wide_data %>%
  filter(tissue == 'buccal',
         pearls_2 == 'high' | pearls_5 == 'high')

# create boostrapping function that outputs matrix of coefficient estimates
################################################################################
boot_sample <- function(i, no_df, high_df, obs_fit) {
  
  # no_df = blood_no
  # high_df = blood_high
  # obs_fit = obs_blood_fit
  
  # take sample from no PEARLS participants
  boot_no <- no_df %>%
    # randomly sample clusters (individuals) with replacement
    slice_sample(n = 10, replace = TRUE) %>%
    # assign new IDs to individuals (since old IDs may be repeated)
    mutate(new_id = 1:10) %>%
    # make data longer, creating separate rows for each timepoint per individual
    pivot_longer(
      cols = !c(pearls_id, new_id, tissue),
      names_to = c(".value", "timepoint"),
      values_to = "timepoint",
      names_pattern = "(.*)_(.*)"
    )

  # take sample from high PEARLS participants (same as above)
  boot_high <- high_df %>%
    slice_sample(n = 10, replace = TRUE) %>%
    # make new IDs different from those for no PEARLS participants
    mutate(new_id = 11:20) %>%
    pivot_longer(
      cols = !c(pearls_id, new_id, tissue),
      names_to = c(".value", "timepoint"),
      values_to = "timepoint",
      names_pattern = "(.*)_(.*)"
    )

  # combine samples from no PEARLS and high PEARLS groups
  boot_all <- bind_rows(boot_no, boot_high)
  # boot_all = boot_sample_df # bootstrap sample that throws error (not warning)
  # obs_fit = obs_blood_fit
  
  lm_gee <- tryCatch(
    expr = {
      warn_msg <- 'none'  # collector for warning message
      
      result <- withCallingHandlers(
        gee(
          horvath2 ~ age*pearls + pearls + age + sex + income_FPL_100,
          id      = new_id,
          data    = boot_all,
          family  = "gaussian",
          corstr  = "exchangeable",
          silent  = TRUE
        ),
        warning = function(w) {
          warn_msg <<- conditionMessage(w)  # save warning, don't stop
          invokeRestart("muffleWarning")    # suppress console output
        }
      )
      
      list(result = result, warning = warn_msg, error = 'none')
    },
    error = function(e) {
      # save problem sample to global env for inspection
      #boot_sample_df <<- boot_all 
      list(result = setNames(as.list(rep(NA, length(coef(obs_fit)))), 
                             names(coef(obs_fit))
                             ), 
           warning = 'none', 
           error = conditionMessage(e))
    }
  )

  # save coefficient estimates and any error/warnings in dataframe
  # extract list of coefficient estimates (or NA if error was generated)
  # coef_estimates <- ifelse(lm_gee$error == 'none',
  #                     as.list(coef(lm_gee$result)),
  #                     as.list(lm_gee$result)
  #                     )
  # 
  # coef_estimates %>%
  #   # convert to dataframe for easier downstream use
  #   data.frame() %>% 
  #   # match variable names to observed model fit
  #   setNames(., names(coef(obs_fit))) %>% 
  #   # add in any error or warnings from model fit
  #   mutate(error = lm_gee$error,
  #          warning = lm_gee$warning
  #   )
  
  if(lm_gee$error == 'none'){
    as.list(coef(lm_gee$result)) %>% # create named list of coefficients
      # convert to dataframe for easier downstream use
      data.frame() %>%
      # match variable names to observed model fit
      setNames(., names(coef(obs_fit))) %>%
      # add in any error or warnings from model fit
      mutate(error = lm_gee$error,
             warning = lm_gee$warning
      )

  } else{
  as.list(lm_gee$result) %>% # create named list of coefficients
    # convert to dataframe for easier downstream use
    data.frame() %>%
    # match variable names to observed model fit
    setNames(., names(coef(obs_fit))) %>%
    # add in any error or warnings from model fit
    mutate(error = lm_gee$error,
           warning = lm_gee$warning
    )
  }
  
}

# set seed for reproducibiility
set.seed(123)

# specify number of bootstrap samples (do 10k at least)
nb <- 10000

# take bootstrap samples and save estimates + errors/warnings in dataframe
# blood samples
blood_boot <- do.call(what = rbind, 
                      args = lapply(X = seq_len(nb), # return list of length nb
                                    FUN = boot_sample, # function
                                    blood_no, # function argument 1
                                    blood_high, # function argument 2
                                    obs_blood_fit # function argument 3
                                    )
                      )

# buccal samples
buccal_boot <- do.call(what = rbind, 
                      args = lapply(X = seq_len(nb), # return list of length nb
                                    FUN = boot_sample, # function
                                    buccal_no, # function argument 1
                                    buccal_high, # function argument 2
                                    obs_buccal_fit
                                    )
                      )

# save matrices
save(blood_boot, buccal_boot, file = "data-processed/bootstrap-estimates.RData")