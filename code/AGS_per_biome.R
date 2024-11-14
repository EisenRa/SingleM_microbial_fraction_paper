library(tidyverse)
library(tinytable)

types_of_interest = c('human gut metagenome','marine metagenome','soil metagenome', 
                      'food metagenome','human oral metagenome','plant metagenome')

df <- read_delim("data/sra/per_acc_summary.csv.gz") %>%
  filter(metagenome_size > 500000000) %>%
  filter(organism %in% types_of_interest) %>%
  mutate(AGS_mb = average_bacterial_archaeal_genome_size / 1000000)

df %>%
  summarise(mean_AGS = mean(AGS_mb),
            median_AGS = median(AGS_mb),
            max_AGS = max(AGS_mb),
            min_AGS = min(AGS_mb),
            .by = organism) %>% 
  tt(., digits = 3)
