# From Ben
library(ggbeeswarm)

c_data <- read_delim("data/pipe/simulated/2pc_read_fractions.tsv")

options(repr.plot.width=6, repr.plot.height=4)

fig2C <- ggplot(data = c_data,  
                aes(x = fct_relevel(novelty_category2, "novel phylum", "novel class", 
                                    "novel order", "novel family", "novel genus", "novel species"),
                    y = (bacterial_archaeal_bases/metagenome_size) * 100)) +
  geom_beeswarm(cex=2) +
  ylab('SingleM microbial fraction (%)') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    title = element_text(size = 15, face = "bold"),
  ) +
  geom_hline(yintercept=(c(1.20,0.8))*100, colour='red') +
  geom_hline(yintercept=100*(c(1.1,0.9)), colour='blue') +
  scale_y_log10() +
  xlab('') +
  ylab("SingleM Microbial Fraction (%)")

ggsave(plot = fig2C, filename = "figures/SI_figureXX.png")
