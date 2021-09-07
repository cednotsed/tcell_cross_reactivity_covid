  ## Visualise unexposed TCR data  ##
  rm(list=ls())
  setwd('../Desktop/git_repos/covid_cross_immunity/')
  require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggplot2); require(ggpubr); require(see)
  
  df <- fread("data/deconvoluted_epitopes/prop_of_unexplained_epitopes.merged.csv")
  
  # prepare text annotations
  annot_df <- df %>%
    filter(Accession %in% c("MN996532", "AY274119", 
                            "MG772933", "MG772934", "EPIISL1098866")) %>%
    mutate(text = case_when(
      Accession == "MN996532" ~ "RaTG13",
      Accession == "MT040335" ~ "PCoV_GX",
      Accession == "AY274119" ~ "SARS Tor2",
      Accession == "EPIISL412977" ~ "RmYN02",
      Accession == "EPIISL1098866" ~ "PrC31",
      Accession == "EPIISL1699443" ~ "RsYN03",
      Accession == "EPIISL1699444" ~ "RsYN04",
      Accession == "EPIISL1699445" ~ "RmYN05",
      Accession == "EPIISL1699446" ~ "RpYN06",
      Accession == "EPIISL1699447" ~ "RmYN07",
      Accession == "EPIISL1699448" ~ "RmYN08",
      Accession == "EPIISL1699449" ~ "RsYN09",
      Accession == "EPIISL852604" ~ "RShSTT182",
      Accession == "EPIISL852605" ~ "RShSTT200",
      Accession == "MG772933" ~ "bat-SL-CoVZC45",
      Accession == "MG772934" ~ "bat-SL-CoVZXC21")) %>%
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
                position = position_nudge(x = -0.1, y = 0),
                angle = 30,
                hjust = 1,
                color = "black", 
                size = 2) +
      labs(x = "Viral genus", y = "Proportion of epitopes explained") +
      coord_flip() +
      theme(legend.position = "none")
  plt
  
  ggsave("results/deconvoluted_epitopes/proportion_of_unexplained_epitopes_explained.merged.png", plot = plt, dpi = 300, width = 8, height = 4)
    
