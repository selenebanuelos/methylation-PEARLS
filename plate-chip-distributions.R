### Author: Selene Banuelos
### Date: 3/09/2026
### Description: Visualize distribution of DNAm age and sample characteristics
### across sample processing plate and methylation array chips
### Note: all samples processed in 1 plate

# setup
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer) # color palettes

# import data
################################################################################
# DNA methylation age estimates
DNAm_age <- read.csv('data-processed/dnam-age-sample-info.csv')

# blood sample information
bl_info <- readRDS('data-raw/Final_SampleInfo_Blood_n39.rds')

# buccal sample information
bu_info <- readRDS('data-raw/Final_SampleInfo_Buccal_n38.rds')

# data wrangling
################################################################################
# specify batch variables
batch_vars <- c('Sample_Well', # sample processing well on plate
                'Sentrix_ID', # chip array ID
                'Sentrix_Position' # chip array row
                )

# clean up batch information
batch <- mutate(bu_info, Sentrix_ID = as.factor(Sentrix_ID)) %>%
  bind_rows(., bl_info) %>% # combine buccal and blood sample info
  select(specimenid, Tissue, ACE, all_of(batch_vars)) %>% # vars to keep
  # separate sample plate column and row information from well ID
  mutate(Sample_Row = substr(Sample_Well, start = 1, stop = 1),
         Sample_Col = as.integer(
           substr(Sample_Well, start = 2, stop = 3)
         )
         ) %>%
  # modify vars for joining downstream
  mutate(specimenid = as.numeric(as.character(specimenid)),
         tissue = stringr::str_to_lower(Tissue)
         ) %>%
  select(-Tissue) # remove duplicate

# combine DNAm age and batch variables
all_data <- full_join(DNAm_age, batch, by = c('specimenid', 'tissue'))

# make sample processing plate data longer for plotting
plate <- all_data %>%
  select(horvath2, 
         ped_be, 
         Sample_Col, 
         Sample_Row, 
         tissue, 
         ACE,
         pearls_id,
         timepoint
         ) %>%
  pivot_longer(cols = c('horvath2', 'ped_be'),
               names_to = 'clock',
               values_to = 'DNAm_age'
               )

chip <- all_data %>%
  select(horvath2, 
         ped_be, 
         Sentrix_ID, 
         Sentrix_Position, 
         tissue, 
         ACE, 
         pearls_id,
         timepoint
         ) %>%
  pivot_longer(cols = c('horvath2', 'ped_be'),
               names_to = 'clock',
               values_to = 'DNAm_age'
  ) %>%
  mutate(row = substr(Sentrix_Position, start = 3, stop = 3),
         chip_id = substr(Sentrix_ID, start = 11, stop = 12),
         PEARLS = case_when(ACE == 'Low' ~ 'No',
                            ACE == 'High' ~ 'High')
         )

# data visualization
################################################################################
# save chips at:
# save plates at:
# specify positions of borders between array chips to add to plots
chip_border_positions <- seq(from = 0.5, to = 10.5, by = 1)

# heatmap to visualize DNAm age distribution across sample processing plate
plate_age_map <- ggplot(plate, aes(Sample_Col, Sample_Row)) +
  geom_tile(aes(fill =  DNAm_age)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  facet_wrap(~clock) +
  labs(title = 'DNAm Age by sample processing plate position',
       x = 'Column',
       y = 'Row'
       ) +
  theme_minimal()
plate_age_map

# heatmap to visualize DNAm age distribution across array chips and chip rows
chip_age_map <- ggplot(chip, aes(chip_id, row)) +
  geom_tile(aes(fill =  DNAm_age)) +
  facet_wrap(~clock) +
  labs(title = 'DNAm Age by chip position',
       x = 'Chip',
       y = 'Row'
  ) +
  geom_vline(xintercept = chip_border_positions, colour = 'white', linewidth = 1.75) +
  theme_minimal() 
chip_age_map

# visualize tissue distribution across sample processing plate
plate_tissue_map <- ggplot(plate, aes(Sample_Col, Sample_Row)) +
  geom_tile(aes(fill =  tissue)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  labs(title = 'Tissue type by sample processing plate position',
       x = 'Column',
       y = 'Row'
  ) +
  theme_minimal()
plate_tissue_map

# visualize tissue distribution across array chips and chip rows
chip_tissue_map <- ggplot(chip, aes(chip_id, row)) +
  geom_tile(aes(fill =  tissue)) +
  labs(title = 'Tissue type by chip position',
       x = 'Chip',
       y = 'Row'
  ) +
  geom_vline(xintercept = chip_border_positions, colour = 'white', linewidth = 1.75) +
  theme_minimal()
chip_tissue_map

# visualize PEARLS status across sample processing plate
plate_PEARLS_map <- ggplot(plate, aes(Sample_Col, Sample_Row)) +
  geom_tile(aes(fill =  ACE)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  labs(title = 'PEARLS status by sample processing plate position',
       x = 'Column',
       y = 'Row'
  ) +
  theme_minimal()
plate_PEARLS_map

# visualize PEARLS status across array chip
chip_PEARLS_map <- ggplot(chip, aes(chip_id, row)) +
  geom_tile(aes(fill = PEARLS)) +
  labs(title = 'PEARLS status by chip position (yellow = no, pink = high',
       x = 'Chip',
       y = 'Row'
  ) +
  scale_fill_manual(breaks = levels(chip$PEARLS),
                    values = c('#ED7673', '#FDF6B5')) +
  geom_vline(xintercept = chip_border_positions, colour = 'white', linewidth = 1.75) +
  theme_minimal()
chip_PEARLS_map

# visualize timepoint distribution across sample processing plate
plate_timepoint_map <- ggplot(plate, aes(Sample_Col, Sample_Row)) +
  geom_tile(aes(fill =  as.factor(timepoint))) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  labs(title = 'Timepoint distribution across sample processing plate',
       x = 'Column',
       y = 'Row'
  ) +
  theme_minimal()
plate_timepoint_map

# visualize timepoint distribution across array chip
chip_timepoint_map <- ggplot(chip, aes(chip_id, row)) +
  geom_tile(aes(fill = as.factor(timepoint))) +
  labs(title = 'Timepoint distribution across chips (light = T2, dark = T5)',
       x = 'Chip',
       y = 'Row'
  ) +
  scale_fill_manual(breaks = levels(chip$PEARLS),
                    values = c('#D4F3A3', '#006D6F')) +
  geom_vline(xintercept = chip_border_positions, colour = 'white', linewidth = 1.75) +
  theme_minimal()
chip_timepoint_map

# visualized subject distribution across chip arrays (all subjects run on one plate)
chip_id_map <- ggplot(chip, aes(chip_id, row)) +
  geom_tile(aes(fill =  pearls_id)) +
  # use dynamic/high contrast colors to visualize many subjects
  scale_fill_manual(values = paletteer_dynamic("cartography::multi.pal", 20)) +
  labs(title = 'Participant distribution on chips',
       x = 'Chip',
       y = 'Row'
  ) +
  geom_vline(xintercept = chip_border_positions, colour = 'white', linewidth = 1.75) +
  theme_minimal()
chip_id_map 