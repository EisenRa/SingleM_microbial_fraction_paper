################################################################################
### Supplementary figure 4
### Raphael Eisenhofer 2024
################################################################################
# import SMF estimates
library(tidyverse)
library(janitor)
library(patchwork)
library(gghighlight)

# Curated metadata from https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001792#sec012
curated_metadata <- read_delim("data/sra/curated_metadata.tsv.gz")
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


#Combine with curated marine metagenome metadata
marine_curated <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  inner_join(., mmdb, by = c("acc" = "library_id")) %>%
  clean_names() %>%
  filter(gbp > 0.5) %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(!is.na(ncbi_stat))


# get mapping rates
#% mapping to catalogue vs singlem estimate for samples from Nishimura & Yoshizawa 2022
#https://www.nature.com/articles/s41597-022-01392-5 [table S1]
nishi_mapping <- readxl::read_xlsx("data/Nishimura_2022_SI_1.xlsx", skip = 2)

prok <- nishi_mapping %>%
  filter(fraction == "prok enriched")

nishi_merged <- marine_curated %>%
  inner_join(., prok, by = join_by(sample_id == sra_sample)) %>%
  select(sample_id, acc, sra_run, singlem_percent, percent_mapped_on_UGCMP, percent_mapped_on_OceanDNA_MAGs, ncbi_stat) %>%
  mutate(percent_mapped_on_OceanDNA_MAGs = as.numeric(percent_mapped_on_OceanDNA_MAGs),
         percent_mapped_on_UGCMP = as.numeric(percent_mapped_on_UGCMP)) %>%
  rename(UGCMP = percent_mapped_on_UGCMP,
         OceanDNA = percent_mapped_on_OceanDNA_MAGs,
         SMF = singlem_percent,
         `NCBI STAT` = ncbi_stat)

nishi_merged %>%
  pivot_longer(., cols = c(SMF, `NCBI STAT`, UGCMP, OceanDNA),
               values_to = "percent") %>%
  ggplot(aes(x = fct_relevel(name, "SMF", "UGCMP", "OceanDNA", "STAT"), 
             y = percent, 
             group = name, 
             colour = name)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.3, alpha = 0.2) +
  theme_classic() + 
  theme(legend.position = 0,
        axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20)) +
  ylab("Percentage (%)")


# calculate DAMR for nishi et al
nishi_merged %>%
  summarise(DAMR = mean(UGCMP / SMF), DAMR_sd = sd(UGCMP / SMF))
