library(tidyverse)

cm2 <- read_delim("data/tera_soil_assembly_checkm2.tsv") %>%
  filter(Completeness > 50) %>%
  mutate(genome_size_corrected = (Genome_Size * (((100 - Completeness) + 100) / 100)) / 1000000)

cm3 <- read_delim("data/tera_soil_assembly_checkm2.tsv") %>%
  filter(Completeness > 50) %>%
  mutate(genome_size_corrected = (Genome_Size / ((Completeness/100) * (1 + (Contamination/100)))) / 1000000)

cm2 %>%
  summarise(mean_AGS = mean(genome_size_corrected),
            median_AGS = median(genome_size_corrected),
            max_GS = max(genome_size_corrected),
            min_GS = min(genome_size_corrected)) %>%
  tinytable::tt(digits = 3)

cm3 %>%
  summarise(mean_AGS = mean(genome_size_corrected),
            median_AGS = median(genome_size_corrected),
            max_GS = max(genome_size_corrected),
            min_GS = min(genome_size_corrected)) %>%
  tinytable::tt(digits = 3)


df %>%
  ggplot(aes(x = read_fraction, y = species_coverage)) +
  geom_point(alpha = 0.3, size = 0.3)

df %>%
  filter(organism != "human gut metagenome" & organism != "gut metagenome") %>%
  ggplot(aes(x = read_fraction, y = species_coverage)) +
  geom_point(alpha = 0.3, size = 0.3) +
  scale_x_continuous(limits = c(0, 100)) + 
  scale_y_continuous(limits = c(0, 1))

df %>%
  filter(organism != "human gut metagenome" & organism != "gut metagenome") %>%
  summarise(n = n(), .by = organism) %>%
  arrange(desc(n))


df_low.spec_high.smf <- df %>%
  filter(read_fraction > 66) %>%
  filter(species_coverage < 0.25)

df_low.spec_high.smf %>%
  ggplot(aes(x = read_fraction, y = species_coverage)) +
  geom_point(alpha = 0.3, size = 0.3) +
  scale_x_continuous(limits = c(0, 100)) + 
  scale_y_continuous(limits = c(0, 1))


df_low.spec_high.smf %>%
  summarise(n = n(), .by = organism) %>%
  arrange(desc(n))


df_soil %>%
  ggplot(aes(x = read_fraction, y = species_coverage)) +
  geom_point(alpha = 0.3, size = 0.3) +
  scale_x_continuous(limits = c(0, 100)) + 
  scale_y_continuous(limits = c(0, 1))


top10 <- df_soil %>% arrange(desc(metagenome_size)) %>%
  slice(1:10) %>%
  select(sample)
