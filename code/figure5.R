################################################################################
### R code for reproducing figure 5.
### Raphael Eisenhofer 03/04/2024
################################################################################


## Load packages + data
library(tidyverse)
library(janitor)
library(patchwork)
library(glue)

acc_organism <- read_delim("data/sra/acc_organism.csv.gz")
ncbi_method <- read_delim("data/sra/NCBI_method_taxonomy_processed.csv.gz")
meta_data_extra <- read_delim("data/sra/extra_metadata_short.tsv.gz") %>%
  clean_names()
singlem_method <- read_delim("data/sra/sra126_r214_mach2.condense.read_fraction6.in_sandpiper.csv.gz",
                             col_names = c("acc", "estimated_bases", "total_bases", "singlem_percent", "warning")) %>%
  mutate(singlem_percent = as.numeric(sub("%", "", singlem_percent))) %>%
  inner_join(., meta_data_extra, by = join_by(acc == run))

#Filter for WGS/RANDOM selection, samples > 0.5 Gbp
merged <- acc_organism %>%
  inner_join(., ncbi_method, by = "acc") %>%
  inner_join(., singlem_method, by = "acc") %>%
  clean_names() %>%
  mutate(ncbi_stat = microbial_fraction * 100) %>%
  filter(library_strategy == "WGS" & library_selection == "RANDOM") %>%
  filter(gbp > 0.5) 


## Plotting
types_of_interest = c('human gut metagenome','marine metagenome','soil metagenome', 
                      'food metagenome','human oral metagenome','plant metagenome')


#Use glue to insert n() for each type
human <- merged %>% 
  filter(organism == "human gut metagenome") %>%
  nrow()
marine <- merged %>% 
  filter(organism == "marine metagenome") %>% 
  nrow()
soil <- merged %>% 
  filter(organism == "soil metagenome") %>% 
  nrow()
food <- merged %>% 
  filter(organism == "food metagenome") %>% 
  nrow()
oral <- merged %>% 
  filter(organism == "human oral metagenome") %>% 
  nrow()
plant <- merged %>% 
  filter(organism == "plant metagenome") %>% 
  nrow()


names_n <- c(paste0('A) human gut metagenome', " (n=", {human}, ")"), paste0('B) human oral metagenome', " (n=", {oral}, ")"),
             paste0('C) marine metagenome', " (n=", {marine}, ")"), paste0('D) soil metagenome', " (n=", {soil}, ")"),
             paste0('E) food metagenome', " (n=", {food}, ")"), paste0('F) plant metagenome', " (n=", {plant}, ")"))
names(names_n) <- c('human gut metagenome', 'human oral metagenome', 'marine metagenome', 
                    'soil metagenome', 'food metagenome', 'plant metagenome')

# stats
# create function to loop through vector
summarise_smf <- function(x){
  merged %>%
  filter(organism == paste0(x)) %>%
  summarise(mean_smf = mean(singlem_percent),
            median_smf = median(singlem_percent),
            sd_smf = sd(singlem_percent)
            )
}

# use map() to loop through vector of names
types_of_interest %>%
  map(., summarise_smf) %>%
  bind_rows %>%
  mutate(metagenome_type = types_of_interest)



# figure 5
merged %>%
  filter(organism %in% types_of_interest) %>%
  mutate(organism = factor(organism, levels = names(names_n))) %>%
  ggplot(aes(x = singlem_percent)) +
  geom_density(size = 2) +
  facet_wrap(~organism, 
             scales = "free_y", 
             labeller = labeller(organism = names_n)) +
  theme_minimal() +
  xlim(0,100) +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 20, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1)
  ) +
  xlab("\nSingleM microbial fraction (%)") +
  ylab("Density")

ggsave("figures/Figure_5.png", width = 15, height = 15, unit = "in")
