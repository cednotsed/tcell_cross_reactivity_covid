## Correlation analysis for expected proportions of hit explained given 
## cophenetic distance from SARS-CoV-2
rm(list=ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggplot2); require(ggpubr); require(see)

tree <- read.tree("data/cov_n2531_ml_tree_delta_rerooted.tree")
prop_df <- fread("data/deconvoluted_epitopes/prop_of_all_epitopes.merged.csv")

cophen_mat <- data.frame(cophenetic.phylo(tree))
cophen_mat$Accession <- rownames(cophen_mat)
cophen_mat <- cophen_mat[, colnames(cophen_mat) %in% c("Accession", "MN908947")]
cophen_mat <- tibble(cophen_mat) %>%
  rename(cophenetic_distance = "MN908947")

plot_df <- cophen_mat %>%
  inner_join(prop_df, by = "Accession") %>%
  mutate(log_cophenetic_distance = log(cophenetic_distance))

# Get unique hosts and species
filtered_df <- plot_df %>%
  distinct(Host, Species_clean, .keep_all = T)

# Linear regression
l1 <- lm(filtered_df$prop_hits ~ filtered_df$log_cophenetic_distance)
filtered_df$studres <- MASS::studres(l1)
corr <- cor.test(y = filtered_df$log_cophenetic_distance, x = filtered_df$prop_hits)

# Annotations
annot_hcov <- filtered_df %>%
  filter(grepl("NL63|OC43|HKU1|229E", Species_clean) & Host == "Homo sapiens")
annot_sars <- filtered_df %>%
  filter(grepl("SARS", Species_clean) & Host == "Homo sapiens")
annot_pcov <- filtered_df %>%
  filter(grepl("pcov", GenBank_Title, ignore.case = T))
annot_ratg13 <- filtered_df %>%
  filter(grepl("ratg13", GenBank_Title, ignore.case = T))
annot_RmYN02 <- filtered_df %>%
  filter(Accession == "EPIISL412977")

corr_p <- signif(corr$p.value, digits = 2)
corr_stat <- round(corr$estimate, 3)

# Plot linear regression + legend
plt1 <- ggplot(filtered_df) +
  geom_point(alpha = 0.3, aes(x = log_cophenetic_distance, y = prop_hits, color = Genus_clean)) +
  geom_smooth(aes(x = log_cophenetic_distance, y = prop_hits), 
              method = "lm", 
              se = T,
              color = "black",
              size = 0.5) + 
  annotate("text", x = 1, y = 0.9,
           label = paste0("R", " = ", corr_stat), 
           hjust = 0,
           size = 3) + 
  annotate("text", x = 1, y = 0.8,
           label = paste0("p", " = ", corr_p),
           hjust = 0,
           size = 3) +
  geom_point(data = annot_hcov, aes(x = log_cophenetic_distance, y = prop_hits),
             color = "black") +
  geom_point(data = annot_sars, aes(x = log_cophenetic_distance, y = prop_hits),
             color = "red") +
  geom_point(data = annot_pcov, aes(x = log_cophenetic_distance, y = prop_hits),
             color = "blue") +
  geom_point(data = annot_ratg13, aes(x = log_cophenetic_distance, y = prop_hits),
             color = "firebrick") +
  geom_point(data = annot_RmYN02, aes(x = log_cophenetic_distance, y = prop_hits),
             color = "saddlebrown") +
  labs(x = "log(distance)", y = "Proportion of epitopes", color = "Viral genus") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Plot linear regression without legend
plt3 <- plt2 +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# Studentised residuals
stud_plt <- ggplot(filtered_df, aes(x = log_cophenetic_distance, y = studres, color = Genus_clean)) +
  geom_point() +
  geom_hline(yintercept = c(-3, 3), color = "red", linetype = "dashed") +
  geom_point(data = annot_hcov, aes(x = log_cophenetic_distance, y = studres),
             color = "black") +
  geom_point(data = annot_sars, aes(x = log_cophenetic_distance, y = studres),
             color = "red") +
  geom_point(data = annot_pcov, aes(x = log_cophenetic_distance, y = studres),
             color = "blue") +
  geom_point(data = annot_ratg13, aes(x = log_cophenetic_distance, y = studres),
             color = "firebrick") +
  geom_point(data = annot_RmYN02, aes(x = log_cophenetic_distance, y = studres),
             color = "saddlebrown") +
  annotate("text", x = c(-1, -1),
            y = c(3.4, -3.3), 
            label = c("y = 3", "y = -3"),
            color = c("red", "red"),
            size = c(3,3)) + 
  labs(x = "ln(distance)", y = "Studentised residuals") +
  theme_bw()

combined <- ggarrange(plt1, stud_plt, nrow = 2, common.legend = T, align = "hv")
combined

ggsave("results/deconvoluted_epitopes/cophenetic_distance_prop.combined.png", plot = combined, dpi = 300, height = 5, width = 8)
fwrite(filtered_df, "results/deconvoluted_epitopes/cophenetic_distance_prop.csv")





