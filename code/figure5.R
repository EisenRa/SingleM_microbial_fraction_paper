################################################################################
## Script to reproduce figure 5. 
## Raphael Eisenhofer Sept 2023
################################################################################

## load libraries/functions
library(tidyverse)
library(patchwork)

read_fraction_files <- list.files("data/pipe/real_data/hyena/r214/",
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

fig5a <- estimates %>%
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

fig5b <- estimates %>%
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

combined_hyenas <- fig5a | fig5b + plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
  )

fig5 <- combined_hyenas + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("figures/Figure_5.png", fig5, width = 15, height = 7.5, unit = "in")
