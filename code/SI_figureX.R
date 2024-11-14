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
singlem_method <- read_delim("data/sra/per_acc_summary.csv.gz") %>%
  rename(acc = sample) %>%
  inner_join(., meta_data_extra, by = join_by(acc == run)) %>%
  mutate(bact_bp_avg = root_coverage * 3400832.56,
         smf_avg = (bact_bp_avg / metagenome_size) * 100)

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



merged_curated %>%
  filter(organism_x == "human gut metagenome") %>%
  summarise(mean_smf = mean(read_fraction),
            mean_tax_free = mean(smf_avg))

homo %>%
  ggplot(aes(x = smf_avg, y = read_fraction)) +
  geom_point()

merged_curated %>%
  filter(organism_x == "soil metagenome") %>%
  summarise(mean_smf = mean(read_fraction),
            mean_tax_free = mean(smf_avg))

aamd_curated %>%
  summarise(mean_smf = mean(read_fraction),
            mean_tax_free = mean(smf_avg))

marine_curated %>%
  summarise(mean_smf = mean(read_fraction),
            mean_tax_free = mean(smf_avg))

