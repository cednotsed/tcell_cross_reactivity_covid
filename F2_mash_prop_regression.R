## Correlation analysis for expected proportions of hit explained given 
## cophenetic distance from SARS-CoV-2
rm(list=ls())
setwd('../Desktop/git_repos/tcell_cross_reactivity_covid')
require(tidyverse);require(data.table);require(ggplot2); require(ggpubr); require(foreach)

prop_df <- fread("data/deconvoluted_epitopes/prop_of_all_epitopes.merged.csv")
mash_df <- distinct(fread("data/mash_distances_k11-21.csv"))

# Linear Regression
linreg <- function(filtered_df) {
  l1 <- lm(filtered_df$prop_hits ~ filtered_df$log_mash_distance)
  filtered_df$studres <- MASS::studres(l1)
  corr <- cor.test(y = filtered_df$log_mash_distance, x = filtered_df$prop_hits)
  corr_p <- signif(corr$p.value, digits = 2)
  corr_stat <- round(corr$estimate, 3)
  
  return(list(df = filtered_df, 
              corr_p = corr_p, 
              corr_stat = corr_stat))
}

# Visualise results
analyse <- function(filtered_df, corr_stat, corr_p) {
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
  
  # Plot linear regression + legend
  plt1 <- ggplot(filtered_df) +
    geom_point(alpha = 0.3, aes(x = log_mash_distance, y = prop_hits, color = Genus_clean)) +
    geom_smooth(aes(x = log_mash_distance, y = prop_hits), 
                method = "lm", 
                se = T,
                color = "black",
                size = 0.5) + 
    annotate("text", x = -2, y = 0.9,
             label = paste0("R", " = ", corr_stat), 
             hjust = 0,
             size = 3) + 
    annotate("text", x = -2, y = 0.8,
             label = paste0("p", " = ", corr_p),
             hjust = 0,
             size = 3) +
    geom_point(data = annot_hcov, aes(x = log_mash_distance, y = prop_hits),
               color = "black") +
    geom_point(data = annot_sars, aes(x = log_mash_distance, y = prop_hits),
               color = "red") +
    geom_point(data = annot_pcov, aes(x = log_mash_distance, y = prop_hits),
               color = "blue") +
    geom_point(data = annot_ratg13, aes(x = log_mash_distance, y = prop_hits),
               color = "firebrick") +
    geom_point(data = annot_RmYN02, aes(x = log_mash_distance, y = prop_hits),
               color = "saddlebrown") +
    labs(x = "log(distance)", y = "Proportion of epitopes", color = "Viral genus") +
    ylim(0, 1) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  # Studentised residuals
  stud_plt <- ggplot(filtered_df, aes(x = log_mash_distance, y = studres, color = Genus_clean)) +
    geom_point() +
    geom_hline(yintercept = c(-3, 3), color = "red", linetype = "dashed") +
    geom_point(data = annot_hcov, aes(x = log_mash_distance, y = studres),
               color = "black") +
    geom_point(data = annot_sars, aes(x = log_mash_distance, y = studres),
               color = "red") +
    geom_point(data = annot_pcov, aes(x = log_mash_distance, y = studres),
               color = "blue") +
    geom_point(data = annot_ratg13, aes(x = log_mash_distance, y = studres),
               color = "firebrick") +
    geom_point(data = annot_RmYN02, aes(x = log_mash_distance, y = studres),
               color = "saddlebrown") +
    annotate("text", x = c(-1, -1),
             y = c(3.4, -3.3), 
             label = c("y = 3", "y = -3"),
             color = c("red", "red"),
             size = c(3,3)) +
    labs(x = "ln(distance)", y = "Studentised residuals") +
    theme_bw()
  
  combined <- ggarrange(plt1, stud_plt, 
                        nrow = 2, 
                        common.legend = T, legend = "bottom",
                        align = "hv")
  
  return(combined)
}
  
# Studentised residuals distribution
plot_studres_dist <- function(filtered_df) {
  studres_plt <- filtered_df %>%
    ggplot(aes(x = studres)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = c(-3, 3), color = "red", linetype = "dashed") +
    annotate("text", y = c(125, 125),
             x = c(3.4, -3.4), 
             label = c("y = 3", "y = -3"),
             color = c("red", "red"),
             size = c(3,3)) +
    labs(y = "Frequency", x = "Studentised residuals")
  return(studres_plt)
}

# MAIN Analysis
foreach(k = seq(11,21), .packages = c("ggplot2")) %do% {
  plot_df <- mash_df %>%
    right_join(prop_df, by = "Accession") %>%
    mutate(log_mash_distance = log(get(str_glue("mash.dist.to.wuhan1.k.{k}"))))
    
  # Get unique hosts and species
  host_species_cophen <- plot_df %>%
    distinct(Host, Species_clean, get(str_glue("mash.dist.to.wuhan1.k.{k}")), .keep_all = T)
  
  # Add Studentised residuals
  HSC_res <- linreg(host_species_cophen)
  HSC_df <- HSC_res[["df"]]
  HSC_corr_stat <- HSC_res[["corr_stat"]]
  HSC_corr_p <- HSC_res[["corr_p"]]
  
  HSC_plt <- analyse(HSC_df, HSC_corr_stat, HSC_corr_p)
  HSC_studres_plt <- plot_studres_dist(HSC_df)
  
  HSC_plt
  HSC_studres_plt
  
  fwrite(HSC_df, 
         str_glue("results/deconvoluted_epitopes/mash_plots/cophenetic_distance_prop.host_species_cophen.k{k}.250821.csv"))
  
  ggsave(str_glue("results/deconvoluted_epitopes/mash_plots/mash_linreg.k{k}.250821.png"),
         plot = HSC_plt,
         dpi = 300,
         height = 5,
         width = 8)
  
  ggsave(str_glue("results/deconvoluted_epitopes/mash_plots/mash_studres.k{k}.250821.png"),
         plot = HSC_studres_plt,
         dpi = 300,
         height = 3,
         width = 8)
  
  nrow(HSC_df)
}

### Visualise table of unexplained epitopes
# ep_list <- read.csv("data/deconvoluted_epitopes/list_of_unexplained_epitopes.txt", header = F)
# 
# tibble(ep_list) %>%
#   separate("V1", into = c("study", "protein", NA)) %>%
#   count(protein)
# 
# ep_df <- tibble(ep_list) %>%
#   separate("V1", into = c("study", "protein", NA))
# 
# ep_tab <- table(ep_df$protein)
# 
# test_df <- tibble(protein = names(ep_tab), n = as.vector(ep_tab)) %>%
#   mutate(core = ifelse(protein == "NSP12"|
#                          protein == "NSP13" |
#                          protein == "NSP14" |
#                          protein == "NSP15" |
#                          protein == "NSP16" |
#                          protein == "N" |
#                          protein == "S", "core", "non-core"))
# test_df %>%
#   ggplot(aes(x = protein, y =  n, fill = core)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Protein", y = "No. of unexplained epitopes", fill = "") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggsave("results/deconvoluted_epitopes/freq_table_unexplained_epitopes.png", width = 5)
