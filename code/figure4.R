################################################################################
## Script to reproduce figure 4. 
## Raphael Eisenhofer Sept 2023
################################################################################

## load libraries/functions
library(tidyverse)
library(patchwork)


read_fraction_files <- list.files("data/pipe/real_data/humans/r214",
                                  pattern = "*_read_fraction.tsv",
                                  full.names = TRUE
                                  )

MAG_mapping_files <- list.files("data/pipe/real_data/humans/",
                                pattern = "*mapping_rate.txt",
                                full.names = TRUE
                                )

assembly_mapping_files <- list.files("data/pipe/real_data/humans/",
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
    catalogue = case_when(
      str_detect(sample, "HRGM") ~ "HRGM",
      !str_detect(sample, "HRGM") ~ "MAGs"
    ),
    sample = str_replace_all(sample, "HRGM_", "")
  ) %>%
  rename("percent_unmapped_MAGs" = percent_unmapped)

estimates <- map_dfr(read_fraction_files, read_delim) %>%
  mutate(sample = str_replace_all(sample, "_M_1", "")) %>%
  left_join(., MAG_mapping_rates, by = "sample") %>%
  left_join(., assembly_mapping_rates, by = "sample") %>%
  mutate(
    dataset = case_when(
      str_detect(catalogue, "HRGM") ~ "HRGM",
      str_detect(sample, "05gb") ~ "human_faecal_0.5gbp",
      str_detect(sample, "_5gb") ~ "human_faecal_5gbp",
      ),
    read_fraction_full = (bacterial_archaeal_bases / metagenome_size) * 100)
    )

#Plotting
colours_humans <- c("human_faecal_0.5gbp" = "#6082B6", "human_faecal_5gbp" = "darkblue", "HRGM" = "orange")

fig4a <- estimates %>%
  ggplot(aes(x = assembly_mapping_percent, y = read_fraction_full, colour = dataset)) +
  geom_point(size = 5, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 102)) +
  scale_color_manual(values = colours_humans) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
  ) +
  labs(y = "SingleM microbial fraction (%)", x = "Percent mapping to assembly (%)") +
  ggtitle("Assembly")

fig4b <- estimates %>%
  ggplot(aes(x = percent_mapped_MAGs, y = read_fraction_full, colour = dataset)) +
  geom_point(size = 5, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 102)) +
  scale_color_manual(values = colours_humans) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
  ) +
  labs(y = "SingleM microbial fraction (%)", x = "Percent mapping to MAG catalogue (%)") +
  ggtitle("MAG catalogue")

combined <- fig4a | fig4b +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
  )

fig4 <- combined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("figures/Figure_4.png", fig4, width = 15, height = 7.5, unit = "in")
