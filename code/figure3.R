################################################################################
## Script to reproduce figure 3. 
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

#Stats for paper
##average % mapping to coassemblies
estimates %>%
  summarise(assembly_mapping = mean(assembly_mapping_percent),
            std_assembly = sd(assembly_mapping_percent),
            mag_mapping = mean(percent_mapped_MAGs),
            std_mag = sd(percent_mapped_MAGs),
            smf_mean = mean(read_fraction_full),
            smf_sd = sd(read_fraction_full),
            .by = dataset)


#Plotting
colours_humans <- c(human_faecal_5gbp = "darkblue", human_faecal_0.5gbp = "grey", HRGM = "orange")


fig3a <- estimates %>%
  ggplot(aes(x = assembly_mapping_percent, y = read_fraction_full, colour = dataset)) +
  geom_point(size = 5, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 102)) +
  scale_color_manual(values = colours_humans) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
  ) +
  labs(y = "SingleM microbial fraction (%)", x = "Percent mapping to assembly (%)") +
  ggtitle("Assembly")

fig3b <- estimates %>%
  ggplot(aes(x = percent_mapped_MAGs, y = read_fraction_full, colour = dataset)) +
  geom_point(size = 5, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 102)) +
  scale_color_manual(values = colours_humans) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
  ) +
  labs(y = "SingleM microbial fraction (%)", x = "Percent mapping to de novo MAGs (%)") +
  ggtitle("de novo MAGs")

combined <- fig3a | fig3b +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    title = element_text(size = 18, face = "bold")
  )

fig3 <- combined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("figures/Figure_3.png", fig4, width = 15, height = 7.5, unit = "in")
