## Visualise unexposed TCR data  ##
rm(list=ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggplot2); require(ggpubr); require(see)

df <- fread("data/deconvoluted_epitopes/prop_of_unexplained_epitopes.merged.csv")

# prepare text annotations
annot_df <- df %>%
  filter(Accession %in% c("MN996532", "MT040335", "AY274119", "EPIISL412977")) %>%
  mutate(text = case_when(
    Accession == "MN996532" ~ "RaTG13",
    Accession == "MT040335" ~ "PCoV_GX",
    Accession == "AY274119" ~ "SARS Tor2",
    Accession == "EPIISL412977" ~ "RmYN02")) %>%
  select(Accession, prop_hits, Genus_clean, text)
  


plt <- ggplot(df, aes(x = Genus_clean, y = prop_hits, color = Genus_clean, fill = Genus_clean)) +
    theme_bw() +
    geom_point(position = position_jitter(width = 0.08), 
               size = 0.5, 
               alpha = 0.5) +
    geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
                 width = 0.2, 
                 outlier.shape = NA,
                 alpha = 0.3) +
    geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
    geom_point(data = annot_df, aes(x = Genus_clean, y = prop_hits), color = "red") +
    geom_text(data = annot_df, aes(x = Genus_clean, y = prop_hits, label = text), 
              position = position_nudge(x = -0.2, y = 0),
              angle = 30,
              hjust = 1,
              color = "black", 
              size = 2) +
    labs(x = "Viral genus", y = "Proportion of epitopes explained") +
    coord_flip() +
    theme(legend.position = "none")
plt

ggsave("results/deconvoluted_epitopes/proportion_of_unexplained_epitopes_explained.merged.png", plot = plt, dpi = 300, width = 8, height = 4)
  
