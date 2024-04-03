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

human_estimates <- map_dfr(read_fraction_files, read_delim) %>%
  mutate(sample = str_replace_all(sample, "_M_1", "")) %>%
  left_join(., MAG_mapping_rates, by = "sample") %>%
  left_join(., assembly_mapping_rates, by = "sample") %>%
  mutate(
    dataset = case_when(
      str_detect(catalogue, "HRGM") ~ "HRGM",
      str_detect(sample, "05gb") ~ "human_faecal_0.5gbp",
      str_detect(sample, "_5gb") ~ "human_faecal_5gbp",
      ),
    read_fraction_full = (bacterial_archaeal_bases / metagenome_size) * 100,
    read_fraction = as.numeric(sub("%", "", read_fraction))
    )

#Stats for paper
##humans
human_estimates %>%
  summarise(assembly_mapping = mean(assembly_mapping_percent),
            std_assembly = sd(assembly_mapping_percent),
            mag_mapping = mean(percent_mapped_MAGs),
            std_mag = sd(percent_mapped_MAGs),
            smf_mean = mean(read_fraction_full),
            smf_sd = sd(read_fraction_full),
            .by = dataset) %>%
  mutate(DAMR = mag_mapping / smf_mean)


# Hyenas
hyena_md <- read_delim("data/msystems.00965-22-s0008.txt") %>%
  unite("sample", hyenaID:sampleID, sep = "_")

read_fraction_files <- list.files("data/pipe/real_data/hyena/r214/",
                                  pattern = "*_read_fraction.tsv.gz",
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

hyena_estimates <- map_dfr(read_fraction_files, read_delim) %>%
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
    read_fraction_full = bacterial_archaeal_bases / metagenome_size * 100,
    read_fraction = as.numeric(sub("%", "", read_fraction)),
    DAMR = percent_mapped_MAGs / read_fraction
  )

#Stats for paper
##hyena
hyena_estimates %>%
  summarise(assembly_mapping = mean(assembly_mapping_percent),
            std_assembly = sd(assembly_mapping_percent),
            mag_mapping = mean(percent_mapped_MAGs),
            std_mag = sd(percent_mapped_MAGs),
            smf_mean = mean(read_fraction),
            smf_sd = sd(read_fraction)) %>%
  mutate(DAMR = mag_mapping / smf_mean)


#Plotting
colours_humans <- c(`de novo MAGs` = "darkblue", human_faecal_0.5gbp = "#488cfa", HRGM = "orange")
colours_hyenas <- c("hyena_D1" = "#ffffcc", "hyena_D3" = "#a1dab4", "hyena_G1" = "#41b6c4", "hyena_G3" = "#225ea8")


fig3a <- human_estimates %>%
  filter(dataset != "human_faecal_0.5gbp") %>%
  mutate(dataset = case_when(str_detect(dataset, "human_faecal_5gbp") ~ "de novo MAGs", .default = "HRGM")) %>%
  ggplot(aes(x = percent_mapped_MAGs, y = read_fraction, colour = dataset)) +
  geom_point(size = 5, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = colours_humans, labels = c("human_faecal_0.5gbp" = "0.5 GBp",
                                                         "human_faecal_5gbp" = "5 GBp")) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    title = element_text(size = 18, face = "bold")
  ) +
  labs(y = "SingleM microbial fraction (%)", x = "Percent mapping to MAGs (%)") +
  ggtitle("Human faecal samples")

fig3b <- hyena_estimates %>%
  ggplot(aes(x = percent_mapped_MAGs, y = read_fraction, fill = dataset)) +
  geom_point(size = 5, alpha = 0.8, shape = 21, stroke = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = colours_hyenas) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    title = element_text(size = 18, face = "bold")
  ) +
  labs(y = "", x = "Percent mapping to MAGs (%)") +
  ggtitle("Hyena faecal samples")

combined <- fig3a | fig3b +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
  )

fig3 <- combined + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("figures/Figure_3.png", fig3, width = 15, height = 7.5, unit = "in")
