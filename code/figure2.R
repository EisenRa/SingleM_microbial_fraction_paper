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
sim_zymo_spiked_homo2_bp_total <- 471442350 * 2
sim_zymo_spiked_arabidopsis_bp_total <- 486523650 * 2
sim_zymo_spiked_arabidopsis2_bp_total <- 365945400 * 2
sim_zymo_spiked_plasmodium_bp_total <- 245668650 * 2
sim_zymo_spiked_plasmodium2_bp_total <- 303076800 * 2

## Total human bp in spiked simulated zymo metagenome?
sim_zymo_bp <- (187352700 * 2)
sim_zymo2_bp <- (186442350 * 2)
sim_zymo_mbp <- sim_zymo_bp/1000000
sim_zymo2_mbp <- sim_zymo2_bp/1000000

## How many bp of bacterial/arcahaeal DNA are there in our Zymo simulated metagenome?
## This is the ground truth, need to minus the fungi and plasmid seqeunces
ecoli_plasmid <- 110007 * 9
salm_plasmid <- 49572 * 9
staph_plasmid <- 11546 * 15
plasmids <- ecoli_plasmid + salm_plasmid + staph_plasmid

ecoli_plasmid2 <- 110007 * 9.1
salm_plasmid2 <- 49572 * 8.9
staph_plasmid2 <- 11546 * 14.9
plasmids2 <- ecoli_plasmid2 + salm_plasmid2 + staph_plasmid2

#N.b. this is # reads in R1 of metagenome * read length * 2
cryptococcus <- 93192 * 150 * 2
saccharomyces <- 42785 * 150 * 2
fungi_bp <- cryptococcus + saccharomyces
fungi2_bp <- (12709200 + 7064100) * 2 #calculated directly from the fungal R1 files

bacterial_archaeal_ground_truth_bp <- sim_zymo_bp - fungi_bp - plasmids
bacterial_archaeal_ground_truth2_bp <- sim_zymo2_bp - fungi2_bp - plasmids2


## Import data
zymo <- read_delim("data/pipe/simulated/Simulated_Zymo_read_fraction.tsv")
zymo_homo <- read_delim("data/pipe/simulated/Simulated_Zymo_Spiked_Homo_read_fraction.tsv")
zymo_plasmodium <- read_delim("data/pipe/simulated/Simulated_Zymo_Spiked_Plasmodium_read_fraction.tsv")
zymo_arabidopsis <- read_delim("data/pipe/simulated/Simulated_Zymo_Spiked_Arabidopsis_read_fraction.tsv")

zymo2 <- read_delim("data/pipe/simulated/Simulated_Zymo2_read_fraction.tsv")
zymo2_homo <- read_delim("data/pipe/simulated/Simulated_Zymo2_Spiked_Homo_read_fraction.tsv")
zymo2_plasmodium <- read_delim("data/pipe/simulated/Simulated_Zymo2_Spiked_Plasmodium_read_fraction.tsv")
zymo2_arabidopsis <- read_delim("data/pipe/simulated/Simulated_Zymo2_Spiked_Arabidopsis_read_fraction.tsv")


## plot (option 2)
zymo_smf <- (zymo$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp) * 100
zymo_homo_smf <- (zymo_homo$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp) * 100
zymo_plasmodium_smf <- (zymo_plasmodium$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp) * 100
zymo_arabidopsis_smf <- (zymo_arabidopsis$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp) * 100

zymo2_smf <- (zymo2$bacterial_archaeal_bases / bacterial_archaeal_ground_truth2_bp) * 100
zymo2_homo_smf <- (zymo2_homo$bacterial_archaeal_bases / bacterial_archaeal_ground_truth2_bp) * 100
zymo2_plasmodium_smf <- (zymo2_plasmodium$bacterial_archaeal_bases / bacterial_archaeal_ground_truth2_bp) * 100
zymo2_arabidopsis_smf <- (zymo2_arabidopsis$bacterial_archaeal_bases / bacterial_archaeal_ground_truth2_bp) * 100

zymo_smf <- (zymo$bacterial_archaeal_bases / zymo$metagenome_size) * 100
zymo_homo_smf <- (zymo_homo$bacterial_archaeal_bases / zymo_homo$metagenome_size) * 100
zymo_plasmodium_smf <- (zymo_plasmodium$bacterial_archaeal_bases / zymo_plasmodium$metagenome_size) * 100
zymo_arabidopsis_smf <- (zymo_arabidopsis$bacterial_archaeal_bases / zymo_arabidopsis$metagenome_size) * 100

zymo2_smf <- (zymo2$bacterial_archaeal_bases / zymo2$metagenome_size) * 100
zymo2_homo_smf <- (zymo2_homo$bacterial_archaeal_bases / zymo2_homo$metagenome_size) * 100
zymo2_plasmodium_smf <- (zymo2_plasmodium$bacterial_archaeal_bases / zymo2_plasmodium$metagenome_size) * 100
zymo2_arabidopsis_smf <- (zymo2_arabidopsis$bacterial_archaeal_bases / zymo2_arabidopsis$metagenome_size) * 100



figure2A_col <- c("Zymo mock" = "#00b0f0", "+ human DNA" = "#e2f0d9", 
                  "+ arabidopsis DNA" = "#51a145", "+ plasmodium DNA" = "#00bfc4")

fig2Aa1 <- as.data.frame(zymo_smf) %>%
  ggplot(aes(x = "Zymo mock 1", 
             y = zymo_smf,
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
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size = 15),
  ) +
  ylab("SingleM microbial fraction (%)") +
  ggtitle("A)") 

fig2Aa2 <- as.data.frame(zymo2_smf) %>%
  ggplot(aes(x = "Zymo mock 2", 
             y = zymo2_smf,
             fill = "#00b0f0")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#00b0f0") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth2_bp / sim_zymo2_bp * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  )



fig2Ab1 <- as.data.frame(zymo_homo_smf) %>%
  ggplot(aes(x = "\\+ *Homo* <br> reads 1", 
             y = zymo_homo_smf,
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
    axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  )

fig2Ab2 <- as.data.frame(zymo2_homo_smf) %>%
  ggplot(aes(x = "\\+ *Homo* <br> reads 2", 
             y = zymo2_homo_smf,
             fill = "#e2f0d9")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#e2f0d9") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth2_bp / sim_zymo_spiked_homo2_bp_total * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  ) 


fig2Ac1 <- as.data.frame(zymo_plasmodium_smf) %>%
  ggplot(aes(x = "\\+ *Plasmodium* <br> reads 1", 
             y = zymo_plasmodium_smf,
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
    axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  ) 

fig2Ac2 <- as.data.frame(zymo2_plasmodium_smf) %>%
  ggplot(aes(x = "\\+ *Plasmodium* <br> reads 2", 
             y = zymo2_plasmodium_smf,
             fill = "#00bfc4")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#00bfc4") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth2_bp / sim_zymo_spiked_plasmodium2_bp_total * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  ) 


fig2Ad1 <- as.data.frame(zymo_arabidopsis_smf) %>%
  ggplot(aes(x = "\\+ *Arabidopsis* <br> reads 1", 
             y = zymo_arabidopsis_smf,
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
    axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  )

fig2Ad2 <- as.data.frame(zymo2_arabidopsis_smf) %>%
  ggplot(aes(x = "\\+ *Arabidopsis* <br> reads 2", 
             y = zymo2_arabidopsis_smf,
             fill = "#51a145")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#51a145") +
  geom_hline(yintercept = bacterial_archaeal_ground_truth2_bp / sim_zymo_spiked_arabidopsis2_bp_total * 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_blank()
  )


fig2A <- fig2Aa1 + fig2Aa2 + fig2Ab1 + fig2Ab2 + fig2Ac1 + fig2Ac2 + fig2Ad1 + fig2Ad2 +
  plot_layout(nrow = 1)


### CAMI2 simulated metagenomes (figure 2B)  

## Import data
cami_marine <- read_delim("data/pipe/simulated/CAMI2_marine0_read_fraction.tsv")
cami_marine1 <- read_delim("data/pipe/simulated/CAMI2_marine1_read_fraction.tsv")
cami_strain <- read_delim("data/pipe/simulated/CAMI2_strain0_read_fraction.tsv")
cami_strain1 <- read_delim("data/pipe/simulated/CAMI2_strain1_read_fraction.tsv")

# work out smf (estimate / ground truth)
marine0_groundtruth_bp <- 2497105500 * 2
marine0_groundtruth_plasmids_bp <- 56371500 * 2
marine0_groundtruth_bacteria_bp <- marine0_groundtruth_bp - marine0_groundtruth_plasmids_bp
smadness0_groundtruth_bp <- 998253450 * 2
marine_0_smf <- (cami_marine$bacterial_archaeal_bases / marine0_groundtruth_bacteria_bp) * 100
smadness_0_smf <- (cami_strain$bacterial_archaeal_bases / smadness0_groundtruth_bp) * 100
cami_marine$bacterial_archaeal_bases / marine0_groundtruth_bp
marine0_groundtruth_bacteria_bp / marine0_groundtruth_bp

marine1_groundtruth_bp <- 2497597650 * 2
marine1_groundtruth_plasmids_bp <- 37587000 * 2
marine1_groundtruth_bacteria_bp <- marine1_groundtruth_bp - marine1_groundtruth_plasmids_bp
smadness1_groundtruth_bp <- 998287650 * 2
marine_1_smf <- (cami_marine1$bacterial_archaeal_bases / marine1_groundtruth_bacteria_bp) * 100
smadness_1_smf <- (cami_strain1$bacterial_archaeal_bases / smadness1_groundtruth_bp) * 100
cami_marine1$bacterial_archaeal_bases / marine1_groundtruth_bp
marine1_groundtruth_bacteria_bp / marine1_groundtruth_bp

cami <- tibble(
  "marine0" = marine_0_smf,
  "strain0" = smadness_0_smf,
  "marine1" = marine_1_smf,
  "strain1" = smadness_1_smf
) %>%
  pivot_longer(cols = everything(), 
               values_to = "percentage", 
               names_to = "metagenome")

#Plot it

fig2Bi <- cami %>%
  filter(str_detect(metagenome, "marine")) %>%
  ggplot(aes(x = metagenome, 
             y = percentage,
             fill = metagenome)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 97.7) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_fill_manual(values = c("darkblue", "darkblue")) +
  theme_classic() + 
  theme(
    legend.position = "none",
    title = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size = 15),
  ) +
  ggtitle("B)")

fig2Bii <- cami %>%
  filter(str_detect(metagenome, "strain")) %>%
  ggplot(aes(x = metagenome, 
             y = percentage,
             fill = metagenome)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 100) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_fill_manual(values = c("black", "black")) +
  theme_classic() + 
  theme(
    legend.position = "none",
    title = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
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
    axis.title.x = element_text(vjust = 14),
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
  plot_layout(widths = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 1, 3))



#fig2 <- fig2top / fig2C

ggsave("figures/Figure_2.png", fig2, width = 10, height = 6, unit = "in")

