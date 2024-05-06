################################################################################
### R code for reproducing figure 4.
### Raphael Eisenhofer 29/9/2023
################################################################################


## Load packages + data
library(tidyverse)
library(janitor)
library(patchwork)
library(gghighlight)

# Curated metadata from https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001792#sec012
curated_metadata <- read_delim("data/sra/curated_metadata.tsv.gz")
aamd <- read_delim("data/sra/hmgdb_selected_dataset_all_col_20231009_135657.csv.gz")
mmdb <- read_delim("data/sra/marmdb_selected_dataset_all_col_20231103_134849.csv.gz")
acc_organism <- read_delim("data/sra/acc_organism.csv.gz")
ncbi_method <- read_delim("data/sra/NCBI_method_taxonomy_processed.csv.gz")
meta_data_extra <- read_delim("data/sra/extra_metadata_short.tsv.gz") %>%
  clean_names()
singlem_method <- read_delim("data/sra/sra126_r214_mach2.condense.read_fraction6.in_sandpiper.csv.gz",
                             col_names = c("acc", "estimated_bases", "total_bases", "singlem_percent", "warning")) %>%
  mutate(singlem_percent = as.numeric(sub("%", "", singlem_percent))) %>%
  inner_join(., meta_data_extra, by = join_by(acc == run))

#Load and filter datasets
merged <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  clean_names() %>%
  mutate(ncbi_stat = microbial_fraction * 100)

#filter only WGS, RANDOM selection, and samples with > 0.5 Gbp
merged_filtered <- merged %>%
  filter(library_strategy == "WGS" & library_selection == "RANDOM") %>%
  filter(gbp > 0.5) %>%
  filter(!is.na(ncbi_stat))

#Combine with human faecal curated metadata
merged_curated <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  inner_join(., curated_metadata, by = c("acc" = "run_accession")) %>%
  clean_names() %>%
  filter(library_strategy == "WGS" & library_selection == "RANDOM") %>%
  filter(gbp > 0.5) %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(!is.na(ncbi_stat))

#Combine with animal metagenome curated metadata
aamd_curated <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  inner_join(., aamd, by = c("acc" = "library_id")) %>%
  clean_names() %>%
  filter(library_strategy == "WGS" & library_selection == "RANDOM") %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(gbp > 0.5) %>%
  filter(!is.na(ncbi_stat))

#Combine with curated marine metagenome metadata
marine_curated <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  inner_join(., mmdb, by = c("acc" = "library_id")) %>%
  clean_names() %>%
  filter(gbp > 0.5) %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(!is.na(ncbi_stat))

# Number of different sample types
merged_filtered %>%
  group_by(organism) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

## number of SingleM warnings
merged_filtered %>%
  summarise(n = n(), .by = warning)


################################################################################
## Plots
################################################################################

# Create common theme
theme_RE <- theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 12),
)

################################################################################
# All data
################################################################################
fig4a <- merged_filtered %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_abline(linewidth = 1) +
  geom_abline(linewidth = 1, intercept = 5, colour = "grey") +
  geom_abline(linewidth = 1, intercept = -5, colour = "grey") +
  geom_point(size = 0.3, alpha = 0.1) +
  theme_minimal() +
  theme_RE +
  labs(y = "SingleM microbial fraction", x = "") +
  ggtitle("All samples")

fig4a <- ggExtra::ggMarginal(fig4a, type="histogram", size=20, fill = "beige")

################################################################################
# Misc STAT vs SingleM stats
################################################################################

## % of time where SingleM estimates are +- 5% of NCBI STAT:
merged_filtered %>%
  filter(singlem_percent - ncbi_stat <= 5 & ncbi_stat - singlem_percent <= 5) %>%
  select(singlem_percent, ncbi_stat) %>%
  nrow() / nrow(merged_filtered)

## % of instances where SingleM estimate > NCBI STAT (by at least 5%):
merged_filtered %>%
  filter(singlem_percent - ncbi_stat >= 5) %>%
  nrow() / nrow(merged_filtered)

## % of instances where singlem 2X> ncbi_stat
merged_filtered %>%
  filter(singlem_percent > (ncbi_stat * 2) ) %>%
  nrow() / nrow(merged_filtered)

## converse, instances where NCBI STAT > SingleM estimates:
merged_filtered %>%
  filter(ncbi_stat > (singlem_percent * 2)) %>%
  nrow() / nrow(merged_filtered)


##odd
ncbi_winners <- merged_filtered %>%
  filter(ncbi_stat > (singlem_percent * 2))

list <- ncbi_winners %>%
  summarise(n = n(), perc = n()/nrow(ncbi_winners), .by = organism)

ncbi_winners %>%
  ggplot(aes(x = singlem_percent)) +
  geom_density() +
  ggtitle("Density of the ~9,000 metagenomes where STAT >2x SingleM")

a<-ncbi_winners %>%
  ggplot(aes(x = ncbi_stat, y = singlem_percent)) +
  geom_point(size = 0.5, alpha = 0.4)

b<-ncbi_winners %>%
  ggplot(aes(x = ncbi_stat)) +
  geom_boxplot() +
  theme_void()

a/b


################################################################################
# Human gut
################################################################################

# N of samples
merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  nrow()

# all stats
merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  filter(!is.na(ncbi_stat)) %>%
  summarise(median = median(singlem_percent),
            mean = mean(singlem_percent),
            sd = sd(singlem_percent),
            median_ncbi = median(ncbi_stat),
            mean_ncbi = mean(ncbi_stat),
            sd_ncbi = sd(ncbi_stat))


#african / south american stats
merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  filter(continent == "Africa" | continent == "South America") %>%
  filter(!is.na(ncbi_stat)) %>%
  summarise(median = median(singlem_percent),
            mean = mean(singlem_percent),
            sd = sd(singlem_percent),
            median_ncbi = median(ncbi_stat),
            mean_ncbi = mean(ncbi_stat),
            sd_ncbi = sd(ncbi_stat),
            .by = continent)

#Mann whitney tests
merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  filter(continent == "Africa") %>%
  filter(!is.na(ncbi_stat)) %>%
  summarise(wilcox = wilcox.test(singlem_percent, ncbi_stat)$p.value)

merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  filter(continent == "South America") %>%
  filter(!is.na(ncbi_stat)) %>%
  summarise(wilcox = wilcox.test(singlem_percent, ncbi_stat)$p.value)


# All human
fig4b <- merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_abline(linewidth = 1) +
  geom_abline(linewidth = 1, intercept = 5, colour = "grey") +
  geom_abline(linewidth = 1, intercept = -5, colour = "grey") +
  geom_point(size = 0.3, alpha = 0.25) +
  theme_minimal() +
  theme_RE +
  theme(axis.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  labs(y = "", x = "") +
  ggtitle("Human gut metagenomes")

fig4b <- ggExtra::ggMarginal(fig4b, type="histogram", size=20, fill = "beige")


# Just african / SAmerican samples
fig4c <- merged_curated %>%
  filter(organism == "human gut metagenome") %>%
  filter(continent == "Africa" | continent == "South America") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat, colour = continent)) +
  geom_abline(linewidth = 1) +
  geom_abline(linewidth = 1, intercept = 5, colour = "grey") +
  geom_abline(linewidth = 1, intercept = -5, colour = "grey") +
  geom_point(size = 1, alpha = 0.5) +
  scale_colour_manual(values = c("Africa" = "blue", "South America" = "darkgreen")) +
  theme_minimal() +
  theme_RE +
  theme(
    legend.position = c(0.7, 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  ) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  labs(y = "", x = "") +
  ggtitle("African/South American gut")

fig4c <- ggExtra::ggMarginal(fig4c, type="histogram", size=20, fill = "beige")


################################################################################
# Animal metagenomes
################################################################################

# aamd mean stats
aamd_curated %>%
  summarise(n = n(), 
            median = median(singlem_percent),
            mean = mean(singlem_percent),
            sd = sd(singlem_percent),
            median_ncbi = median(ncbi_stat),
            mean_ncbi = mean(ncbi_stat),
            sd_ncbi = sd(ncbi_stat)
            )

# aamd stats by taxon
aamd_taxon <- aamd_curated %>%
  summarise(n = n(), 
            median = median(singlem_percent),
            mean = mean(singlem_percent),
            sd = sd(singlem_percent),
            median_ncbi = median(ncbi_stat),
            mean_ncbi = mean(ncbi_stat),
            sd_ncbi = sd(ncbi_stat),
            .by = taxon_name)


#AAMDB
fig4d <- aamd_curated %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_abline(linewidth = 1) +
  geom_abline(linewidth = 1, intercept = 5, colour = "grey") +
  geom_abline(linewidth = 1, intercept = -5, colour = "grey") +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  theme_RE +
  labs(y = "SingleM microbial fraction", x = "STAT microbial fraction") +
  ggtitle("Non-human animal metagenomes")

fig4d <- ggExtra::ggMarginal(fig4d, type="densigram", size=20, fill = "beige")


#stats
aamd_curated %>%
  summarise(mean_singlem = mean(singlem_percent, na.rm = T),
            sd_singlem = sd(singlem_percent, na.rm = T),
            mean_ncbi = mean(ncbi_stat, na.rm = T),
            sd_ncbi = sd(ncbi_stat, na.rm = T))

################################################################################
# Marine metagenome samples
################################################################################

fig4e <- marine_curated %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_abline(linewidth = 1) +
  geom_abline(linewidth = 1, intercept = 5, colour = "grey") +
  geom_abline(linewidth = 1, intercept = -5, colour = "grey") +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  theme_RE +
  labs(y = "", x = "STAT microbial fraction") +
  ggtitle("Marine metagenomes")

fig4e <- ggExtra::ggMarginal(fig4e, type="histogram", size=20, fill = "beige")


################################################################################
# Soil metagenomes
################################################################################
fig4f <- merged_filtered %>%
  filter(organism == "soil metagenome") %>%
  ggplot(aes(y = singlem_percent, x = ncbi_stat)) +
  geom_abline(linewidth = 1) +
  geom_abline(linewidth = 1, intercept = 5, colour = "grey") +
  geom_abline(linewidth = 1, intercept = -5, colour = "grey") +
  geom_point(size = 0.75, alpha = 0.5) +
  theme_minimal() +
  theme_RE +
  theme(axis.title.y = element_blank()) +
  labs(y = "SingleM microbial fraction", x = "STAT microbial fraction") +
  ggtitle("Soil metagenomes")

fig4f <- ggExtra::ggMarginal(fig4f, type="histogram", size=20, fill = "beige")




################################################################################
# Create composite figure 4
################################################################################
Figure_4 <- wrap_elements(fig4a) + wrap_elements(fig4b) + wrap_elements(fig4c) + wrap_elements(fig4d) + wrap_elements(fig4e) + wrap_elements(fig4f) + 
  plot_layout(ncol = 3) + plot_annotation(tag_levels = "A")

ggsave("figures/Figure_4.png", Figure_4, width = 12, height = 9, unit = "in")


