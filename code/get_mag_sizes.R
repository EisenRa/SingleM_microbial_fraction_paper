################################################################################
## Script to get hyena MAG sizes for updating read_fraction estimates
## Also making SI figure X
## Raphael Eisenhofer Oct 2023
################################################################################

## load libraries/functions
library(tidyverse)
library(patchwork)

## load data + filter
taxonomy <- read_delim("data/pipe/real_data/hyena/hyena_gtdbtk.tsv") %>%
  select(user_genome, classification) %>%
  mutate(
    classification = str_replace_all(classification, ";.__$", ""),
    classification = str_replace_all(classification, ";.__$", ""),
    classification = str_replace_all(classification, ";.__$", ""),
    classification = str_replace_all(classification, ";.__$", ""),
    classification = str_replace_all(classification, ";.__$", ""),
    rank = word(classification, -1, sep = ";"),
    user_genome = str_replace_all(user_genome, ".fa", "")
    )

#size = genome size scaled to 100% by checkm completeness (not perfect, but will do)
mag_sizes <- read_delim("data/pipe/real_data/hyena/hyena_metawrap_70_10_bins.stats") %>%
  select(bin, size, completeness) %>%
  mutate(size = size * (((100 - completeness) + 100) / 100))

# merge tables + create summary
combined <- taxonomy %>%
  inner_join(mag_sizes, by = join_by(user_genome == bin)) %>%
  select(rank, size) %>%
  summarise(
    genome_size = mean(size),
    min_genome_size = min(size),
    max_genome_size = max(size),
    .by = rank
    )

write_delim(combined, "data/pipe/real_data/hyena/hyena_mag_genome_sizes.tsv", delim = "\t")


### BASH code
#cut -f1 hyena_mag_genome_sizes.tsv | sed '1d;' > hyena_ranks.tsv
#grep -vwf hyena_ranks.tsv r214_gtdb_mean_genome_sizes.tsv > sans_hyena_ranks.tsv
#sed '1d;' hyena_mag_genome_sizes.tsv > hyena_mag_genome_sizes_no_header.tsv
#cat sans_hyena_ranks.tsv hyena_mag_genome_sizes_no_header.tsv > r214_hyena_mag_sizes.tsv



## Comparison of values when MAG sizes are used
read_fraction_files_magsizes <- list.files("data/pipe/real_data/hyena/r214_magsizes/",
                                  pattern = "*_read_fraction.tsv",
                                  full.names = TRUE
)

read_fraction_files <- list.files("data/pipe/real_data/hyena/r214/",
                                  pattern = "*_read_fraction.tsv",
                                  full.names = TRUE
)

rf1 <- map_dfr(read_fraction_files, read_delim) %>%
  mutate(sample = str_replace_all(sample, "_M_1", ""),
         database = "r214",
         read_fraction_full = bacterial_archaeal_bases / metagenome_size)
rf2 <- map_dfr(read_fraction_files_magsizes, read_delim) %>%
  mutate(sample = str_replace_all(sample, "_M_1", ""),
         database = "r214_mag_sizes",
         read_fraction_full = bacterial_archaeal_bases / metagenome_size)

comparison <- rbind(rf1, rf2)

comparison %>%
  ggplot(aes(x = sample, y = read_fraction_full, colour = database, group = sample)) +
  geom_point(size = 4) +
  geom_line(colour = "black") +
  theme_classic() +
  geom_hline(yintercept = 1.0) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 90)
  )

comparison %>%
  ggplot(aes(x = database, y = read_fraction_full, colour = database)) +
  geom_boxplot() +
  geom_point(size = 4) +
  theme_classic() +
  geom_hline(yintercept = 1.0) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14)
  )

##############
## SI figure X
##############

read_fraction_files <- list.files("data/pipe/real_data/hyena/r214_magsizes/",
                                  pattern = "*_read_fraction.tsv",
                                  full.names = TRUE
)

MAG_mapping_files <- list.files("data/pipe/real_data/hyena/",
                                pattern = "*mapping_rate.txt",
                                full.names = TRUE
)

assembly_mapping_files <- list.files("data/pipe/real_data/hyena/",
                                     pattern = "*summary.tsv",
                                     full.names = TRUE
)

assembly_mapping_rates <- map_dfr(assembly_mapping_files, read_delim) %>%
  select(sample, assembly_mapping_percent) %>%
  mutate(sample = str_replace_all(sample, "_M", ""))


MAG_mapping_rates <- map_dfr(MAG_mapping_files, read_delim) %>%
  filter(Genome == "unmapped") %>%
  pivot_longer(., !Genome,
               names_to = "sample",
               values_to = "percent_unmapped",
               values_drop_na = TRUE
  ) %>%
  select(!Genome) %>%
  mutate(
    sample = str_replace_all(sample, "_M Relative Abundance \\(%\\)", ""),
    percent_mapped_MAGs = 100 - percent_unmapped,
  ) %>%
  rename("percent_unmapped_MAGs" = percent_unmapped)

estimates <- map_dfr(read_fraction_files, read_delim) %>%
  mutate(sample = str_replace_all(sample, "_M_1", "")) %>%
  left_join(., MAG_mapping_rates, by = "sample") %>%
  left_join(., assembly_mapping_rates, by = "sample") %>%
  mutate(
    dataset = case_when(
      str_detect(sample, "D1_") ~ "hyena_D1",
      str_detect(sample, "D3_") ~ "hyena_D3",
      str_detect(sample, "G1_") ~ "hyena_G1",
      str_detect(sample, "G3_") ~ "hyena_G3"
    ),
    read_fraction_full = bacterial_archaeal_bases / metagenome_size * 100
  )

#Plotting
colours_hyenas <- c("hyena_D1" = "#ffffcc", "hyena_D3" = "#a1dab4", "hyena_G1" = "#41b6c4", "hyena_G3" = "#225ea8")

si_fig_Xa <- estimates %>%
  ggplot(aes(x = assembly_mapping_percent, y = read_fraction_full, fill = dataset)) +
  geom_point(size = 5, alpha = 0.8, shape = 21, stroke = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 120)) +
  scale_y_continuous(limits = c(0, 120)) +
  scale_color_manual(values = colours_hyenas) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
  ) +
  labs(y = "SingleM microbial fraction (%)", x = "Percent mapping to assembly (%)") +
  ggtitle("Assemblies")

si_fig_Xb <- estimates %>%
  ggplot(aes(x = percent_mapped_MAGs, y = read_fraction_full, fill = dataset)) +
  geom_point(size = 5, alpha = 0.8, shape = 21, stroke = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 120)) +
  scale_y_continuous(limits = c(0, 120)) +
  scale_color_manual(values = colours_hyenas) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
  ) +
  labs(y = "SingleM microbial fraction (%)", x = "Percent mapping to MAG catalogue (%)") +
  ggtitle("MAG catalogue")

combined_hyenas <- si_fig_Xa | si_fig_Xb + plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
  )

si_fig_X <- combined_hyenas + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("figures/Figure_SX.png", si_fig_X, width = 15, height = 7.5, unit = "in")
