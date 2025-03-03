################################################################################
## Script to reproduce figure 2. 
## Raphael Eisenhofer / Ben Woodcroft Nov 2023
################################################################################

## load libraries/functions
library(tidyverse)
library(patchwork)

### Simulated zymo mock communities with spiked eukaryotic DNA

## Total size of the spiked simulated zymo metagenomes (base pairs)?
sim_zymo_spiked_homo_bp_total <- 757352700 * 2
sim_zymo_spiked_arabidopsis_bp_total <- 486523650 * 2
sim_zymo_spiked_plasmodium_bp_total <- 245668650 * 2

## Total human bp in spiked simulated zymo metagenome?
sim_zymo_bp <- (187352700 * 2)
sim_zymo_mbp <- sim_zymo_bp/1000000

human_bp <- sim_zymo_spiked_homo_bp_total - sim_zymo_bp
arabidopsis_bp <- sim_zymo_spiked_arabidopsis_bp_total - sim_zymo_bp
plasmodium_bp <- sim_zymo_spiked_plasmodium_bp_total - sim_zymo_bp

## How many bp of bacterial/arcahaeal DNA are there in our Zymo simulated metagenome?
## This is the ground truth, need to minus the fungi and plasmid seqeunces
ecoli_plasmid <- 110007 * 9
salm_plasmid <- 49572 * 9
staph_plasmid <- 11546 * 15
plasmids <- ecoli_plasmid + salm_plasmid + staph_plasmid

#N.b. this is # reads in R1 of metagenome * read length * 2
cryptococcus <- 93192 * 150 * 2
saccharomyces <- 42785 * 150 * 2
fungi_bp <- cryptococcus + saccharomyces

bacterial_archaeal_ground_truth_bp <- sim_zymo_bp - fungi_bp - plasmids


## Import data
zymo <- read_delim("data/pipe/simulated/Simulated_Zymo_read_fraction.tsv")
zymo_homo <- read_delim("data/pipe/simulated/Simulated_Zymo_Spiked_Homo_read_fraction.tsv")
zymo_plasmodium <- read_delim("data/pipe/simulated/Simulated_Zymo_Spiked_Plasmodium_read_fraction.tsv")
zymo_arabidopsis <- read_delim("data/pipe/simulated/Simulated_Zymo_Spiked_Arabidopsis_read_fraction.tsv")

## plot (option 2)
zymo_accuracy <- (zymo$bacterial_archaeal_bases / sim_zymo_bp) * 100
zymo_homo_accuracy <- (zymo_homo$bacterial_archaeal_bases / sim_zymo_spiked_homo_bp_total) * 100
zymo_plasmodium_accuracy <- (zymo_plasmodium$bacterial_archaeal_bases / sim_zymo_spiked_plasmodium_bp_total) * 100
zymo_arabidopsis_accuracy <- (zymo_arabidopsis$bacterial_archaeal_bases / sim_zymo_spiked_arabidopsis_bp_total) * 100

figure2A_col <- c("Zymo mock" = "#00b0f0", "+ human DNA" = "#e2f0d9", 
                  "+ arabidopsis DNA" = "#51a145", "+ plasmodium DNA" = "#00bfc4")

fig2Aa <- as.data.frame(zymo_accuracy) %>%
  ggplot(aes(x = "Zymo mock", 
             y = zymo_accuracy,
             fill = "#00b0f0")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#00b0f0") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth_bp / sim_zymo_bp * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size = 15),
  ) +
  ylab("SingleM microbial fraction (%)") +
  ggtitle("A)") 

#accuracy:
zymo_accuracy / 100 / (bacterial_archaeal_ground_truth_bp / sim_zymo_bp) 
zymo_homo_accuracy / 100 / (bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_homo_bp_total) 
zymo_arabidopsis_accuracy / 100 / (bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_arabidopsis_bp_total) 
zymo_plasmodium_accuracy / 100 / (bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_plasmodium_bp_total)

#accuracy (total bact/arch)
zymo$read_fraction
(bacterial_archaeal_ground_truth_bp / sim_zymo_bp) * 100
zymo_homo$read_fraction
(bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_homo_bp_total) * 100
zymo_arabidopsis$read_fraction
(bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_arabidopsis_bp_total) * 100
zymo_plasmodium$read_fraction
(bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_plasmodium_bp_total) * 100


fig2Ab <- as.data.frame(zymo_homo_accuracy) %>%
  ggplot(aes(x = "+ Homo \nreads", 
             y = zymo_homo_accuracy,
             fill = "#e2f0d9")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#e2f0d9") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_homo_bp_total * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  )


fig2Ac <- as.data.frame(zymo_plasmodium_accuracy) %>%
  ggplot(aes(x = "\\+ *Plasmodium* <br> reads", 
             y = zymo_plasmodium_accuracy,
             fill = "#00bfc4")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#00bfc4") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_plasmodium_bp_total * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 15, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  ) 

fig2Ad <- as.data.frame(zymo_arabidopsis_accuracy) %>%
  ggplot(aes(x = "\\+ *Arabidopsis* <br> reads", 
             y = zymo_arabidopsis_accuracy,
             fill = "#51a145")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#51a145") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth_bp / sim_zymo_spiked_arabidopsis_bp_total * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 15, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  )

fig2A <- fig2Aa + fig2Ab + fig2Ac + fig2Ad +
  plot_layout(nrow = 1)


### CAMI2 simulated metagenomes (figure 2B)  

## Import data
cami_marine <- read_delim("data/pipe/simulated/CAMI2_marine0_read_fraction.tsv")
cami_strain <- read_delim("data/pipe/simulated/CAMI2_strain0_read_fraction.tsv")

# work out accuracy (estimate / ground truth)
marine0_groundtruth_bp <- 2497105500 * 2
marine0_groundtruth_plasmids_bp <- 56371500 * 2
marine0_groundtruth_bacteria_bp <- marine0_groundtruth_bp - marine0_groundtruth_plasmids_bp
smadness0_groundtruth_bp <- 998253450 * 2
marine_0_accuracy <- (cami_marine$bacterial_archaeal_bases / marine0_groundtruth_bacteria_bp) * 100
smadness_0_accuracy <- (cami_strain$bacterial_archaeal_bases / smadness0_groundtruth_bp) * 100
cami_marine$bacterial_archaeal_bases / marine0_groundtruth_bp
marine0_groundtruth_bacteria_bp / marine0_groundtruth_bp

cami <- tibble(
  "marine" = marine_0_accuracy,
  "strain" = smadness_0_accuracy
) %>%
  pivot_longer(cols = everything(), 
               values_to = "percentage", 
               names_to = "metagenome")

#Plot it

fig2Bi <- cami %>%
  filter(metagenome == "marine") %>%
  ggplot(aes(x = metagenome, 
             y = percentage,
             fill = metagenome)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 97.7) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_fill_manual(values = c("darkblue")) +
  theme_classic() + 
  theme(
    legend.position = "none",
    title = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size = 15),
  ) +
  ggtitle("B)")

fig2Bii <- cami %>%
  filter(metagenome == "strain") %>%
  ggplot(aes(x = metagenome, 
             y = percentage,
             fill = metagenome)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_fill_manual(values = c("black")) +
  theme_classic() + 
  theme(
    legend.position = "none",
    title = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 

fig2B <- fig2Bi | fig2Bii


#Panel C

#import data
file_paths <- list.files(path = "benchmarks/increasing_novelty_benchmark/output_singlem/singlem/", full.names = T) %>%
  str_sort(., numeric = T)

novelty_df <- read_delim(paste0(file_paths, "/marine0.smf")) %>%
  mutate(novel = seq(0, 100, 10),
         smf_uncapped = (bacterial_archaeal_bases / metagenome_size) * 100)

#plot plot
fig2C <- novelty_df %>%
  ggplot(aes(x = novel, y= smf_uncapped)) + 
  geom_point() + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.title.x = element_text(vjust = 25),
    title = element_text(size = 15, face = "bold"),
  ) +
  ylim(70,120) + 
  geom_line() + 
  labs(x = "% of community\nthat is novel") + 
  geom_hline(yintercept=c(110, 90), colour='blue') + 
  geom_hline(yintercept=c(120, 80), colour='red') +
  geom_hline(yintercept=c(100), colour='grey') +
  ggtitle("C)")



fig2 <- fig2A + fig2B + fig2C +
  plot_layout(widths = c(0.5, 0.5, 0.5, 0.5, 1, 3))

#fig2 <- fig2top / fig2C

ggsave("figures/Figure_2.png", fig2, width = 10, height = 6, unit = "in")

