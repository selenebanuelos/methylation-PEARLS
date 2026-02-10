### Author: Selene Banuelos
### Date: 12/15/2025
### Description: Assess all chronological age predictor performance by 
### calculating correlation between predicted age and chronological age as well
### as mean absolute error between predicted and chronological age

# setup
library(tidyverse)
library(rmcorr) # calculate correlation for repeated measures

# import data
################################################################################
# information on required CpGs per clock
clock_info <- read.csv('data-processed/summary_methscore_CpG.csv')

# generated methylation predictions
predictions <- read.csv('data-processed/methylation-predictions.csv')

# sample information
sample <- read.csv('data-processed/pearls-acesmatchingbysexage.csv')

# data wrangling 
################################################################################
# names of chronological age predictors
chrono_age <- c('horvath_age', 
                'horvath2', 
                'hannum_age', 
                'epigenetic_age_zhang', 
                'c_age', 
                'ped_be'
                )

# calculate % of missing predictive CpGs and filter for chronological age clocks
clocks <- clock_info %>%
  # format predictor name for joining downstream
  mutate(predictor = janitor::make_clean_names(predictor),
         # calculate % of missing CpGs
         perc_missing_cpg = 100*(nCpG_missing / nCpG_required)
  ) %>%
  filter(predictor %in% chrono_age)

# clean up DNAme sample info for joining with methylation predictions 
sample_long <- sample %>%
  pivot_longer(cols = c(T2_specimenid, T5_specimenid),
               names_pattern = '(.)_(.*)',
               names_to = c('timepoint', '.value')
               ) %>%
  mutate(specimenid = as.character(specimenid)) %>%
  dplyr::rename(pearls_id = subjectid)

# standardize format of clock names to match chrono_age vec and for joining
names(predictions) <- janitor::make_clean_names(names(predictions))

# join methylation predictions with sample information
joint <- predictions %>%
  mutate(specimenid = str_remove(sample_id,'^0+' # remove leading zeros
                                 ),
         specimenid = str_remove(specimenid,'..$' # remove last two characters
                                 )
         ) %>%
  # keep predictions, sample id, chorno age, and vars for filtering
  select(chrono_age, specimenid, age, imp_method, tissue) %>%
  # add in participant ids and participant info
  full_join(., 
            sample_long, 
            by = 'specimenid'
            )

# calculate & plot correlations between epigenetic and chronological ages 
################################################################################
# calculate correlation within blood sample data
blood_corr_horvath2 <- joint %>%
  #filter(timepoint == 2) %>%
  filter(tissue == 'blood') %>%
  rmcorr(participant = 'pearls_id',
         measure1 = 'horvath2',
         measure2 = 'age',
         dataset = .)

blood_corr_pedbe <- joint %>%
  #filter(timepoint == 2) %>%
  filter(tissue == 'blood') %>%
  rmcorr(participant = 'pearls_id',
         measure1 = 'ped_be',
         measure2 = 'age',
         dataset = .)
  
# calculate pearson correlation within buccal sample data
buccal_corr_horvath2 <- joint %>%
  #filter(timepoint == 2) %>%
  filter(tissue == 'buccal') %>%
  rmcorr(participant = 'pearls_id',
         measure1 = 'horvath2',
         measure2 = 'age',
         dataset = .)

buccal_corr_pedbe <- joint %>%
  #filter(timepoint == 2) %>%
  filter(tissue == 'buccal') %>%
  rmcorr(participant = 'pearls_id',
         measure1 = 'ped_be',
         measure2 = 'age',
         dataset = .)

# bar plots showing correlation between epigenetic & chrono age in blood
blood_corr %>%
  # format numbers
  mutate(corr_chrono_age = signif(corr_chrono_age, digits = 2)
  ) %>%
  ggplot(aes(
    x = clock, 
    y = corr_chrono_age)) +
  geom_bar(stat = 'identity') +
  # label correlations between ages within bars
  geom_text(aes(label = corr_chrono_age),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
            ) +
  ylab('Correlation') +
  xlab('Clock') +
  ggtitle('Blood samples: Correlations between epigenetic and chronological age')
  
# bar plots showing correlation between epigenetic & chrono age in buccal swabs
buccal_corr %>%
  # format numbers
  mutate(corr_chrono_age = signif(corr_chrono_age, digits = 2)
  ) %>%
  ggplot(aes(
    x = clock, 
    y = corr_chrono_age)) +
  geom_bar(stat = 'identity') +
  # label correlations between ages within bars
  geom_text(aes(label = corr_chrono_age),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Correlation') +
  xlab('Clock') +
  ggtitle('Buccal samples: Correlations between epigenetic and chronological age')

# calculate & plot mean absolute error between epigenetic and chronological ages 
################################################################################
blood_error <- joint %>%
  filter(tissue == 'blood') %>%
  select(
    c(specimenid, 
      age, 
      all_of(chrono_age)
      )
    ) %>%
  # calculate absolute error for each predictor
  mutate_at(
    # do for all variables except identifiers and chrono age
    vars(!c(specimenid, age)),
    # add '_error' suffix to new variables with absolute error
    list(error = ~ abs( . - age))
    ) %>%
  select(contains('error')) %>% # keep absolute error columns
  pivot_longer(
    cols = everything(),
    names_to = 'predictor',
    values_to = 'abs_error'
  ) %>%
  group_by(predictor) %>%
  # calculate mean absolute error per predictor
  summarise(mean_abs_error = mean(abs_error)) %>%
  mutate(predictor = str_remove(predictor,'_error$'))%>%# remove '_error' suffix
  # add in % missing predictive CpGs
  left_join(.,
            select(clocks, c(predictor, perc_missing_cpg)),
            by = 'predictor')

buccal_error <- joint %>%
  filter(tissue == 'buccal') %>%
  select(
    c(specimenid, 
      age, 
      all_of(chrono_age)
    )
  ) %>%
  # calculate absolute error for each predictor
  mutate_at(
    # do for all variables except identifiers and chrono age
    vars(!c(specimenid, age)),
    # add '_error' suffix to new variables with absolute error
    list(error = ~ abs( . - age))
  ) %>%
  select(contains('error')) %>% # keep absolute error columns
  pivot_longer(
    cols = everything(),
    names_to = 'predictor',
    values_to = 'abs_error'
  ) %>%
  group_by(predictor) %>%
  # calculate mean absolute error per predictor
  summarise(mean_abs_error = mean(abs_error)) %>%
  mutate(predictor = str_remove(predictor,'_error$'))%>%# remove '_error' suffix
  # add in % missing predictive CpGs
  left_join(.,
            select(clocks, c(predictor, perc_missing_cpg)),
            by = 'predictor')

# create bar plots showing mean absolute error 
blood_error %>%
  # format numbers
  mutate(mean_abs_error = signif(mean_abs_error, digits = 2)) %>%
  ggplot(aes(
    x = reorder(predictor, perc_missing_cpg), 
    y = mean_abs_error
  )) +
  geom_bar(stat = 'identity') +
  # label correlations between ages within bars
  geom_text(aes(label = mean_abs_error),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Mean absolute error') +
  xlab('Clock') +
  ggtitle('Blood samples: Mean absolute error between epigenetic and chronological age')

buccal_error %>%
  # format numbers
  mutate(mean_abs_error = signif(mean_abs_error, digits = 2)) %>%
  ggplot(aes(
    x = reorder(predictor, perc_missing_cpg), 
    y = mean_abs_error
  )) +
  geom_bar(stat = 'identity') +
  # label correlations between ages within bars
  geom_text(aes(label = mean_abs_error),
            vjust = 1.5,
            position = position_dodge(width = 0.9)
  ) +
  ylab('Mean absolute error') +
  xlab('Clock') +
  ggtitle('Buccal samples: Mean absolute error between epigenetic and chronological age')

# output 
################################################################################
write.csv(joint,
          file = 'data-processed/dnam-age-sample-info.csv',
          row.names = FALSE
          )