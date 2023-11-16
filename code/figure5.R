################################################################################
### R code for reproducing figure 5.
### Raphael Eisenhofer 29/9/2023
################################################################################

# Curated metadata from https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001792#sec012

## Load packages + data
library(tidyverse)
library(janitor)
library(patchwork)
library(gghighlight)

curated_metadata <- read_delim("data/sra/curated_metadata.tsv.gz")
aamd <- read_delim("data/sra/hmgdb_selected_dataset_all_col_20231009_135657.csv.gz")
mmdb <- read_delim("data/sra/marmdb_selected_dataset_all_col_20231103_134849.csv.gz")
acc_organism <- read_delim("data/sra/acc_organism.csv.gz")
ncbi_method <- read_delim("data/sra/NCBI_method_taxonomy_processed.csv.gz")
singlem_method <- read_delim("data/sra/sra126_r214_mach1.condense.read_fraction7.in_sandpiper.csv.gz",
                             col_names = c("acc", "estimated_bases", "total_bases", "singlem_percent", "warning")) %>%
  mutate(singlem_percent = as.numeric(sub("%", "", singlem_percent)))


#filtering samples with singlem_percent > 100, as these are likely not metagenomes (e.g. single cell amplifications)
merged <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  clean_names() %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(singlem_percent < 100)

merged_curated <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  inner_join(., curated_metadata, by = c("acc" = "run_accession")) %>%
  clean_names() %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(singlem_percent < 100)

aamd_curated <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  inner_join(., aamd, by = c("acc" = "library_id")) %>%
  clean_names() %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(singlem_percent < 100)

marine_curated <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  inner_join(., mmdb, by = c("acc" = "library_id")) %>%
  clean_names() %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(singlem_percent < 100)

# Number of different sample types
merged %>%
  group_by(organism) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

merged_curated %>%
  group_by(organism) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

theme_RE <- theme(
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 14),
)


## Plots

# All data
fig5a <- merged %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 0.5, alpha = 0.1) +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme_RE +
  theme() +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("All samples")

merged_curated %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.2) +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme_RE +
  theme() +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("All samples")

## number of instances where SingleM estimate > NCBI stat:
singlem_greater <- merged %>%
  filter(singlem_percent > ncbi_stat) %>%
  nrow()

singlem_greater / nrow(merged) * 100

singlem_greater2 <- merged %>%
  filter(singlem_percent > (ncbi_stat *2) ) %>%
  nrow()

singlem_greater2 / nrow(merged) * 100

# Human gut
merged %>%
  filter(organism == "human gut metagenome") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.3) +
  #  geom_smooth() +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme_RE +
  theme() +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Human gut metagenomes")

merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.3) +
  #  geom_smooth() +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme_RE +
  theme() +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Human gut metagenomes")

# N of samples
merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  nrow()

# Highlighting non-Western samples
fig5b <- merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat, colour = continent)) +
  geom_point(size = 2, alpha = 0.5) +
  gghighlight(continent == "Africa" | continent == "South America",
              label_key = type
  ) +
  scale_colour_manual(values = c("Africa" = "blue", "South America" = "gold")) +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme_RE +
  theme(
    legend.position = c(0.80, 0.65),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Human gut metagenomes")

# Also for country
merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat, colour = country)) +
  geom_point(size = 2, alpha = 0.5) +
  gghighlight(country == "India" | country == "Mongolia",
              label_key = type
  ) +
  scale_colour_manual(values = c("India" = "blue", "Mongolia" = "gold")) +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme_RE +
  theme(
    legend.position = c(0.15, 0.90),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Human gut metagenomes")


# N with greater singlem > ncbi stat
human_n <- merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  filter(singlem_percent > ncbi_stat) %>%
  nrow()

human_n / merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  nrow()


# Now try to get non-human faecal samples

## These bioprojects appear to be human (though labelled as 'gut metagenome')
# PRJEB34871 = pooled human faecal samples
# PRJEB23010 is also derived from humans
# PRJEB23147 = human colon
# PRJNA496479 is also human
# PRJEB25514 also
# PRJEB14847
# PRJNA518912 = isolates
# PRJEB28237 also cultures
# PRJEB29152 too
# PRJEB8094 humans
# PRJNA217052 humans (fiji)
# PRJEB18687 gut virome
# PRJEB31009 humans
# PRJEB25210 humuan poop pits
# PRJNA421881 human
# PRJNA526822 human
# PRJNA320015
# PRJEB20308

actually_human_samples <- c("PRJEB20308", "PRJNA320015", "PRJNA526822", "PRJNA421881", "PRJEB25210", "PRJEB31009", "PRJEB18687", "PRJEB34871", "PRJEB23010", "PRJEB23147", "PRJNA496479", "PRJEB25514", "PRJEB14847", "PRJEB32631", "PRJNA518912", "PRJEB28237", "PRJEB29152", "PRJEB8094", "PRJNA217052")

merged_curated %>%
  filter(organism == "gut metagenome" | organism == "chicken gut metagenome" | organism == "pig gut metagenome" | organism == "mouse gut metagenome" | organism == "bovine gut metagenome") %>%
  filter(host != "Homo sapiens") %>%
  filter(!bioproject %in% actually_human_samples) %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  theme_RE +
  coord_cartesian(expand = FALSE) +
  theme() +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Putatively non-human gut metagenomes")

# N of samples
merged_curated %>%
  filter(organism == "gut metagenome" | organism == "chicken gut metagenome" | organism == "pig gut metagenome" | organism == "mouse gut metagenome" | organism == "bovine gut metagenome") %>%
  filter(host != "Homo sapiens") %>%
  filter(!bioproject %in% actually_human_samples) %>%
  nrow()

#AAMDB
fig5c <- aamd_curated %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  theme_RE +
  coord_cartesian(expand = FALSE) +
  theme() +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Non-human animal metagenomes")



mar <- merged %>%
  filter(organism == "marine metagenome") %>%
  filter(abs(singlem_percent - ncbi_stat) <= 1) %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  theme_RE 

mar %>%
  summarise(n = n(), .by = bioproject) %>%
  arrange(desc(n))

# marine metagenome
merged %>%
  filter(organism == "marine metagenome") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  theme_RE +
  coord_cartesian(expand = FALSE) +
  theme(
    axis.title.y = element_blank()
  ) +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Marine metagenomes")


##MMDB
marine_curated %>%
  summarise(n = n(), .by = mar_mdb_biome) %>%
  arrange(desc(n))

fig5d <- marine_curated %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  theme_RE +
  coord_cartesian(expand = FALSE) +
  theme(
    axis.title.y = element_blank()
  ) +
  labs(y = "SingleM microbial fraction", x = "NCBI STAT microbial fraction") +
  ggtitle("Marine metagenomes")


# N samples
merged %>%
  filter(organism == "marine metagenome") %>%
  nrow()


Figure_6 <- fig5a + fig5b + fig5c + fig5d + plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")

ggsave("figures/Figure_6.png", Figure_6, width = 15, height = 15, unit = "in")


## number of SingleM warnings
merged %>%
  summarise(n = n(), .by = warning)

4822 / nrow(merged)




hi <- merged %>%
  mutate(warning = str_replace_all(warning, "WARNING.*", "WARNING")) %>%
  filter(warning == "WARNING") 

## Explore enrichment of samples that have singlem > ncbi_stat

singlem_greater <- merged_curated %>%
  mutate(singlem_greater = case_when(singlem_percent > ncbi_stat ~ "TRUE",
         .default = "FALSE"))

singlem_50 <- merged_curated %>%
  mutate(singlem_greater = case_when(singlem_percent > 50 & ncbi_stat < 50 ~ "TRUE",
                                     .default = "FALSE"))


z <- singlem_50 %>%
  summarize(TRUE_count = sum(singlem_greater == TRUE, na.rm = TRUE),
            FALSE_count = sum(singlem_greater == FALSE, na.rm = TRUE),
            .by = organism) 

