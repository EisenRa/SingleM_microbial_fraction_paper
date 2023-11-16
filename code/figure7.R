################################################################################
### R code for reproducing figure 7.
### Raphael Eisenhofer 29/9/2023
################################################################################

# Curated metadata from https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001792#sec012

## Load packages + data
library(tidyverse)
library(janitor)
library(patchwork)
library(glue)

curated_metadata <- read_delim("data/sra/curated_metadata.tsv.gz")
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


## Plotting
types_of_interest = c('human gut metagenome','marine metagenome','soil metagenome', 
                      'food metagenome','human oral metagenome','plant metagenome')

#Use glue to insert n() for each type

human <- merged %>% filter(organism == "human gut metagenome") %>% nrow()
marine <- merged %>% filter(organism == "marine metagenome") %>% nrow()
soil <- merged %>% filter(organism == "soil metagenome") %>% nrow()
food <- merged %>% filter(organism == "food metagenome") %>% nrow()
oral <- merged %>% filter(organism == "human oral metagenome") %>% nrow()
plant <- merged %>% filter(organism == "plant metagenome") %>% nrow()


names_n <- c(paste0('human gut metagenome', " (n=", {human}, ")"), paste0('marine metagenome', " (n=", {marine}, ")"),
             paste0('soil metagenome', " (n=", {soil}, ")"), paste0('food metagenome', " (n=", {food}, ")"),
             paste0('human oral metagenome', " (n=", {oral}, ")"), paste0('plant metagenome', " (n=", {plant}, ")"))
names(names_n) <- c('human gut metagenome','marine metagenome','soil metagenome', 
                    'food metagenome','human oral metagenome','plant metagenome')

merged %>%
  filter(organism %in% types_of_interest) %>%
  ggplot(aes(x = singlem_percent)) +
  geom_density(size = 2) +
  facet_wrap(~organism, scales = "free_y", labeller = labeller(organism = names_n)) +
  theme_minimal() +
  xlim(0,100) +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 20, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("\nSingleM microbial fraction (%)") +
  ylab("Density")

ggsave("figures/Figure_7.png", width = 15, height = 15, unit = "in")
