## plot (option 1)
# work out accuracy (estimate / ground truth)
zymo_accuracy <- (zymo$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp * 100) - 100
zymo_homo_accuracy <- (zymo_homo$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp * 100) - 100
zymo_plasmodium_accuracy <- (zymo_plasmodium$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp * 100) - 100
zymo_arabidopsis_accuracy <- (zymo_arabidopsis$bacterial_archaeal_bases / bacterial_archaeal_ground_truth_bp * 100) - 100
acc_combined <- tibble("Zymo mock" = zymo_accuracy,
                       "+ human DNA" = zymo_homo_accuracy,
                       "+ plasmodium DNA" = zymo_plasmodium_accuracy,
                       "+ arabidopsis DNA" = zymo_arabidopsis_accuracy) %>%
  pivot_longer(., cols = everything(), names_to = "metagenome", values_to = "percentage")

figure2A_col <- c("Zymo mock" = "#00b0f0", "+ human DNA" = "#e2f0d9",
                  "+ arabidopsis DNA" = "#51a145", "+ plasmodium DNA" = "#00bfc4")

level_order <- c("Zymo mock", "+ human DNA", "+ plasmodium DNA", "+ arabidopsis DNA")

fig2A <- acc_combined %>%
  ggplot(aes(x = factor(metagenome, level = level_order),
             y = percentage,
             fill = metagenome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = figure2A_col) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    axis.title.x = element_blank()
  ) +
  ylab("percentage from ground truth (%)")