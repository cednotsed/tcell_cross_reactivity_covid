## Visualise unexposed TCR data  ##
rm(list=ls())
setwd('../Desktop/git_repos/tcell_cross_reactivity_covid')
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggplot2); require(ggpubr)

# Function for getting maximum homology from blast output
parse_blastout <- function(blastout) {
  homology_df <- tibble(Query = epitope_order)
  
  for (hcov in c("229E", "OC43", "NL63", "HKU1")) {
    max_homology_df <- blastout %>%
      rename(Query = qseqid) %>%
      left_join(protein_meta, by = "Query") %>%
      filter(Host == "Homo sapiens" & grepl(hcov, Species_clean)) %>%
      mutate(pident = ifelse(is.na(pident), 0, pident)) %>%
      group_by(Query) %>%
      summarise(max = max(pident, na.rm = T))
    
    colnames(max_homology_df)[2] <- paste0(hcov)
    
    homology_df <- homology_df %>%
      left_join(max_homology_df, by = "Query")
  }
  
  hcov_plot_df <- homology_df %>%
    gather(key = hCoV, value = value, -Query)
  
  hcov_plot_df$Query <- factor(hcov_plot_df$Query, levels = epitope_order)
  hcov_plot_df <- hcov_plot_df %>% mutate(value = ifelse(is.na(value), 0, value))
  
  # Get maximum across four HCoVs
  four_hcov_df <- hcov_plot_df %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>%
    group_by(Query) %>%
    summarise(value = max(value, na.rm = T)) %>%
    mutate(hCoV = "All endemic HCoVs")
  
  hcov_plot_df <- bind_rows(hcov_plot_df, four_hcov_df)
  
  return(hcov_plot_df)
}

# Plot maximum homologies of blastoutput
plot_figure <- function(hcov_plot_df, figure_path) {
  # Computing proportion of 'unexplained' epitopes
  computations <- hcov_plot_df %>% 
    group_by(hCoV) %>%
    # summarise(queries = n_distinct(Query)) # Check if all epitopes are present
    summarise(proportion = sum(value == 0) / nrow(protein_meta)) %>%
    mutate(proportion = round(proportion, 3))
  
  # Arrange order of rows
  hcov_plot_df$hCoV <- factor(hcov_plot_df$hCoV, 
                              levels = rev(c("229E", "NL63", "HKU1", "OC43", "All endemic HCoVs")))
  
  hcov_plot <- hcov_plot_df %>%
    mutate(value = ifelse(value == 0, NA, value)) %>%
    ggplot() +
      geom_tile(aes(y = hCoV, fill = value, x = Query), color = "darkgrey", width=0.9, height=0.95,) +
      geom_vline(xintercept = len_df$cum_len, linetype="dashed", color = "black", size = 0.3) +
      theme(axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 1.8),
            panel.grid = element_blank(), 
            panel.background = element_rect(fill = "white"),
            legend.position = "left") +
      labs(fill = "Max. Homology (%)") +
      scale_fill_gradient(low = "darkseagreen1", high = "darkgreen", na.value = "white") +
      geom_text(data = computations, aes(x = 183, y = hCoV, label = proportion), size = 2) +
      coord_cartesian(xlim = c(0, 185), ylim = c(1, 4.9), # This focuses the x-axis on the range of interest
                      clip = 'off')
  
  ggsave(figure_path, plot = hcov_plot, height = 4, width = 8, dpi = 600)
  
  return(hcov_plot)
}

## Get number of protein queries for annotating figure
protein_meta <- fread("data/deconvoluted_epitopes/deconvoluted_epitopes.csv")
protein_meta$Query <- paste(protein_meta$study, protein_meta$protein, protein_meta$start, sep = "_")

# Get coordinates for annotating proteins on heatmap
protein_order <- c("NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP10", 
                   "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", 
                   "S", "ORF3", "E", "M", "ORF6", "ORF7", "ORF8", "N")

epitope_order <- c()
len <- c()

for (p in protein_order) {
  epitopes <- unique(protein_meta$Query[protein_meta$protein == p])
  epitope_order <- c(epitope_order, epitopes)
  len <- c(len, length(epitopes))
}

len_df <- data.frame(protein = protein_order, len = len)
len_df$cum_len <- 0

# Make len_df cumulative
for (i in seq(nrow(len_df))) {
  len_df[i, "cum_len"] <- sum(len_df[1:i, "len"]) + 0.5
}

len_df <- len_df[-nrow(len_df), ]  # Remove last row

## Analyse data ###############################################################
meta <- fread("data/parsed_host_metadata_FB_LvD_CT.n2572.csv")
figure_path <- "results/deconvoluted_epitopes/deconvoluted_hcov_max_homology.merged.png"
blastout <- fread("data/deconvoluted_epitopes/parsed_deconvoluted_blastout.tsv")

blastout <- blastout %>%
  left_join(meta, by = "Accession")

# Calculate max homology from blast output
hcov_plot_df <- parse_blastout(blastout = blastout)

# Visualise max homology of epitopes to HCoVs
plot_figure(hcov_plot_df = hcov_plot_df, figure_path = figure_path)

# Calculate percentage of unexplainable epitopes that can be explained by each virus
get_epitope_prop <- function(epitope_list, path_to_output) {
  result_df <- tibble()
  
  to_test <- meta %>%
    filter(Species_clean != "SARS-CoV-2") %>%
    filter(Accession %in% unique(blastout$Accession))
  
  to_test <- to_test$Accession
  
  for (acc in to_test) {
    hits <- blastout %>%
      select(qseqid, pident, Accession) %>%
      mutate(pident = ifelse(is.na(pident), 0, pident)) %>%
      filter(Accession == acc) %>%
      filter(qseqid %in% epitope_list) %>%
      group_by(qseqid) %>%
      summarise(max = max(pident, na.rm = T))
  
    n_hits <- nrow(hits)
    prop_hits <- n_hits / length(epitope_list)
    result_df <- rbind(result_df, tibble(Accession = c(acc), prop_hits = c(prop_hits)))
  }
  
  final_results <- result_df %>%
    left_join(meta, by = "Accession") %>%
    arrange(desc(prop_hits))
  
  fwrite(final_results, path_to_output)
  return(final_results)
  
}

# Find epitopes 'unexplained' by HCoVs
unexplained <- hcov_plot_df %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  filter(hCoV != "All endemic HCoVs") %>% # Not necessary but for cleanliness
  group_by(Query) %>%
  summarise(max = max(value, na.rm = T)) %>%
  filter(max == 0)

unexplained <- as.vector(unexplained$Query)

unexplained_df <- get_epitope_prop(epitope_list = unexplained, 
                                   path_to_output = "data/deconvoluted_epitopes/prop_of_unexplained_epitopes.merged.csv")

fwrite(tibble(unexplained), 
       "data/deconvoluted_epitopes/list_of_unexplained_epitopes.txt", 
       row.names = F, 
       col.names = F)
# Frequency table of unexplained epitopes
tibble(unexplained) %>% 
  separate(unexplained, into = c("study", NA, NA), sep = "_") %>%
  group_by(study) %>%
  summarise(n = n())