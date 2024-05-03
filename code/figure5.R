################################################################################
### R code for reproducing figure 5.
### Raphael Eisenhofer 03/04/2024
################################################################################


## Load packages + data
library(tidyverse)
library(janitor)
library(patchwork)
library(glue)
library(emo)
#note I needed to set my graphics device backend to 'AGG': tools -> global options -> graphics

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


names_n <- c(paste0('human gut ', emo::ji("poop"), " (n=", {human}, ")"), paste0('human oral ', emo::ji("tooth"), " (n=", {oral}, ")"),
             paste0('marine ', emo::ji("ocean"), " (n=", {marine}, ")"), paste0('soil ', emo::ji("new_moon"), " (n=", {soil}, ")"),
             paste0('food ', emo::ji("bread"), " (n=", {food}, ")"), paste0('plant ', emo::ji("leaf"), " (n=", {plant}, ")"))
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
sll <- read_delim("code/soil/smf_vs_latitude.generated.csv")

tropic = 23.43614

fig5a <- ggplot(sll) +
  geom_point(aes(x=latitude, y=smf, colour=longitude), alpha=0.3) +
  geom_smooth(aes(x=latitude, y=smf), method='loess') +
  coord_cartesian(expand = F) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(0.1, 0.2),
    legend.background = element_rect()
  ) +
  geom_vline(xintercept=c(0, -90, 90)) +
  geom_vline(xintercept=c(-tropic,tropic), linetype='dashed') +
  scale_x_continuous(breaks = c(-80, -40, 0, 40, 80),
                     labels = c("-80°", "-40°", "0°", "40°", "80°")) +
  scale_colour_gradient(limits = c(-180, 180), 
                        breaks = c(-180, -90, 0, 90, 180),
                        labels = c("-180°", "-90°", "0°", "90°", "180°")
                        ) +
  labs(x = "Latitude", y = "SingleM microbial fraction")


fig5b <- merged %>%
  filter(organism %in% types_of_interest) %>%
  mutate(organism = factor(organism, levels = names(names_n))) %>%
  ggplot(aes(x = singlem_percent)) +
  geom_density(size = 2, fill = "beige") +
  facet_wrap(~organism, 
             scales = "free_y", 
             labeller = labeller(organism = names_n)) +
  coord_cartesian(expand = c(0), clip = "off") +
  theme_minimal() +
  xlim(0,100) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = NA, linewidth = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1)
  ) +
  xlab("\nSingleM microbial fraction") +
  ylab("Density")


fig5a / fig5b +
  plot_layout(heights = c(1, 0.66)) +
  plot_annotation(tag_levels = "A")

ggsave("figures/Figure_5.png", width = 12.5, height = 10, unit = "in")
