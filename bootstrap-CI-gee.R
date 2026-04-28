### Author: Selene Banuelos
### Date: 4/5/2025
### Description: Bootstrap 95% confidence interval for regression coefficient
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
# set seed for reproducibiility
set.seed(123)

# specify number of bootstrap samples
nb <- 10000

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

# specify number of columns in bootstrap coef matrix (# of covariates + 1 for error codes)
p <- length(coef(obs_blood_fit)) + 1

# create empty matrix to save coefficients from each bootstrap sample
coef_mat <- matrix(0, nb, p)

# need to name matrix with variable names corresponding to coefficients
colnames(coef_mat) <- c(names(obs_blood_fit$coefficients), 'error_code')

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
boot_sample <- function(i, no_df, high_df) {

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
  print(head(boot_no))
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
  print(head(boot_high))
  # combine samples from no PEARLS and high PEARLS groups
  boot_all <- bind_rows(boot_no, boot_high)
  
  # fit linear model using generalized estimating equations
  lm_gee <- tryCatch(
    gee(
      horvath2 ~ age*pearls + pearls + age + sex + income_FPL_100,
      id = new_id, # make sure to use newly assigned IDs
      data = boot_all,
      family = 'gaussian',
      corstr = 'exchangeable',
      silent = TRUE
      ),
    error = function(e) NULL
  )
  
  # save coefficient of interest and any error codes
  if (is.null(lm_gee)) return(rep(NA_real_, p))
  c(lm_gee$coefficients, lm_gee$error)
  
}

# take nb bootstrap samples and save coefficients in matrix
blood_mat <- coef_mat[,] <- do.call(rbind, lapply(seq_len(nb), boot_sample, blood_no, blood_high))
buccal_mat <- coef_mat[,] <- do.call(rbind, lapply(seq_len(nb), boot_sample, buccal_no, buccal_high))

# save matrices
save(blood_mat, buccal_mat, file = "data-processed/bootstrap-estimates.RData")

# check boostrap distributions for skew
################################################################################
check_bias <- function(boot_mat, # matrix of bootstrapped estimates
                       obs_fit # model object fit with observed data
                       ){

  # calculate bias: mean(bootstrap estimates) - observed estimate
  point <- round(
    obs_fit$coefficients["pearlshigh"], 
    digits = 3
    )
  cat("Point estimate:", point, "\n")
  
  boot_mean <- round(
    mean(boot_mat[, "pearlshigh"], na.rm = TRUE), 
    digits = 3
    )
  cat("Bias (mean bootsrap estimate - observed point estimate):", boot_mean - point, "\n")
  cat("Skewness: sign tells you which tail is longer\n")
  
  # visualize bootstrap coefficient estimate distributions
  boot_mat %>%
    as.data.frame() %>%
    # remove estimates from any resamples that produced errors when fitting model
    filter(V7 == 0) %>%
    select(-V7) %>%
    # make data longer for plotting
    pivot_longer(
      cols = everything(),
      names_to = "Covariate",
      values_to = "Estimate"
    ) %>%
    # plot distributions of all covariates
    ggplot(aes(x = Estimate)) +
    geom_density() +
    facet_wrap(~Covariate)
}

check_bias(blood_mat, obs_blood_fit)
check_bias(buccal_mat, obs_buccal_fit)
# little bit of bias in both bootstrapped distributions of estimate

# calculate bias-corrected and accelerated confidence intervals 
################################################################################
# can't use coxed::bca() for this reason: https://stackoverflow.com/questions/55401615/r-calculate-bca-from-vector-of-bootstrapped-results
# use rms::bootBCa() to generate BCa bootstrap on existing boostrap replicates
library(rms)

# get coefficient estimates from obsrved sample
obs_est_blood <- round(obs_blood_fit$coefficients, digits = 3)
obs_est_buccal <- round(obs_buccal_fit$coefficients, digits = 3)

# create matrix with only bootstrap estimates (remove error column)
boot_est_blood <- blood_mat[, colnames(blood_mat) != '']
boot_est_buccal <- buccal_mat[, colnames(buccal_mat) != '']

# generate BCa 95% CI that corrects for skew in bootstrapped distribution
bootBCa(estimate = obs_est_blood, 
        estimates = boot_est_blood,
        type = 'bca',
        n = 39, # need to double check this, is this asking for # of independent clusters or total observations
        seed = .Random.seed, # uses random seed set before bootstrap. where is this being used?
        conf.int = 0.95
)
print(obs_est_blood)

bootBCa(estimate = obs_est_buccal, 
        estimates = boot_est_buccal,
        type = 'bca',
        n = 39, # need to double check this, is this asking for # of independent clusters or total observations
        seed = .Random.seed, # uses random seed set before bootstrap. where is this being used?
        conf.int = 0.95
)
print(obs_est_buccal)

### I don't think these are correct since it's unclear what to supply for the 
### 'n' argument - total observations in original sample or the number of 
### independent clusters. Asked Alan and he said he'll look into it. He said that
### for now, just use 2.5th and 97.5th quantiles of the bootstrap distribution 
### to get the 95% CI
### also need to double check if i should remove any bootstrap estimates from 
### model fits that generated error codes