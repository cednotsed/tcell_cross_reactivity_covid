rm(list=ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggplot2); require(ggpubr); require(see)

hit_list <- list.files("results/non_CoV_hits/1000_alignments/", pattern = "hit")  
desc_list <- list.files("results/non_CoV_hits/1000_alignments/", pattern = "desc")
merged_df <- tibble()

for (i in seq(length(hit_list))) {
  df <- fread(paste0("results/non_CoV_hits/1000_alignments/", desc_list[i]))
  df <- df %>% 
    separate(Accession, sep = ",", into = c(NA, "sseqid")) %>%
    mutate(sseqid = gsub('\\"|\\)', "", sseqid))
  
  hit <- fread(paste0("results/non_CoV_hits/1000_alignments/", hit_list[i]))
  colnames(hit) <- c("qseqid", "sseqid", "pident", 
                     "length", "mismatches", "gapopen", 
                     "qstart", "qend", "sstart", "send", 
                     "eval", "bitscore", "perc_pos")
  
  print(ncol(hit))
  print(ncol(df))
  
  temp_df <- hit %>%
    left_join(df, by = "sseqid") %>%
    select(qseqid, sseqid, Description, 
           "Scientific Name", "Common Name", pident, 
           "Query Cover", mismatches, "eval", 
           sstart, send, qstart, qend)
  
  merged_df <- bind_rows(merged_df, temp_df)
}

merged_df %>%
  filter(grepl("schulien", qseqid))
merged_df <- merged_df %>%
  filter(!grepl("synthetic|SARS|coronavirus", Description)) %>%
  rename(qcov = "Query Cover") %>%
  mutate(qcov = str_replace_all(qcov, "[[:punct:]]", "")) %>%
  mutate(qcov = as.numeric(qcov))

protein_meta <- fread("data/deconvoluted_epitopes/deconvoluted_epitopes.csv")
protein_meta <- protein_meta %>% 
  unite("qseqid", study, protein, start, sep = "_")

merged_df <- protein_meta %>%
  right_join(merged_df, by = "qseqid") %>%
  arrange(desc(qcov, pident))

merged_df <- merged_df %>%
  rename(sci_name = "Scientific Name") %>%
  mutate(split = strsplit(sci_name, " "), n_words = lengths(split))

parsed_df <- merged_df %>%
  filter(sci_name != "") %>%
  separate(sci_name, into = c("tax_group", "tax_group2"), sep = " ", remove = F) %>%
  mutate(tax_group = ifelse(tax_group == "Human", "Alphapapillomavirus", tax_group),
         tax_group = ifelse(tax_group %in% c("Candidatus", "unclassified", "uncultured"), tax_group2, tax_group)) %>%
  select(-split, -n_words) %>%
  filter(tax_group != "bacterium" & sci_name != "Gallus gallus") %>%  # exclude uncultured bacterium
  filter(qcov > 99)

parsed_df %>%
  # distinct(qseqid)
#   summarise(count = n_distinct(qseqid))
#   View()
  
fwrite(parsed_df, "results/non_CoV_hits/1000_alignments/combined_results.clean.csv")

summarised_df <- parsed_df %>%
  group_by(tax_group) %>%
  summarise(n_queries = n_distinct(qseqid)) %>%
  arrange(desc(n_queries))

fwrite(summarised_df, "results/non_CoV_hits/1000_alignments/summarised_results.clean.csv")
