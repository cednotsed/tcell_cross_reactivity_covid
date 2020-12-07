# Script to parse and merge Batch Entrez + Prokka blast outputs
rm(list=ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(tidyverse);require(data.table);require(ape);require(seqinr);require(adegenet)

## Parse and merge blastout function ###################################################
parse_blastout <- function(path_to_prokka, path_to_entrez, path_to_output) {
  # Prokka
  prokka <- fread(path_to_prokka)
  colnames(prokka) <- c("qseqid","sseqid","pident","length","mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore")
  blastout_NC <- prokka %>%
    filter(grepl("NC", sseqid)) %>%
    separate(sseqid, sep = "\\_", c("acc1", "acc2", NA), extra = "merge", remove = F) %>%
    mutate(Accession = paste0(acc1, "_", acc2)) %>%
    select(-acc1, -acc2)
  
  blastout_rest <- prokka %>%
    filter(!grepl("NC", sseqid)) %>%
    separate(sseqid, sep = "\\_", c("Accession", NA), extra = "merge", remove = F)
  
  prokka <- bind_rows(blastout_NC, blastout_rest)
  
  # Batch Entrez
  entrez <- fread(path_to_entrez)
  colnames(entrez) <- c('qseqid','sseqid','pident','length','mismatch',
                          'gapopen', 'qstart', 'qend', 'sstart', 'send',
                          'evalue', 'bitscore')
  entrez <- entrez %>%
    separate(sseqid, sep = "\\.", c("Accession", NA, NA), extra = "merge", remove = F) %>%
    mutate(qseqid = str_replace(qseqid, "-", "_"))
  
  
  # Merge and get top blast hits
  merged <- bind_rows(prokka, entrez) %>%
    group_by(qseqid, Accession) %>%
    summarise(pident = max(pident, na.rm = T))
  
  fwrite(merged, path_to_output)
  
  return(merged)
}

## For proteome wide ##
# proteome <- parse_blastout(path_to_prokka = "data/wuhan-hu-1/blast_wuhan_epitopes.prokka_db.tsv",
#                            path_to_entrez = "data/wuhan-hu-1/blast_wuhan_epitopes.batch_entrez.tsv",
#                            path_to_output = "data/wuhan-hu-1/parsed_proteome_blastout.tsv")

# ## For 161 deconvoluted epitopes ##
deconv <- parse_blastout(path_to_prokka = "data/deconvoluted_epitopes/blast_deconvoluted_epitopes.prokka_db.tsv",
                         path_to_entrez = "data/deconvoluted_epitopes/blast_deconvoluted_epitopes.batch_entrez.tsv",
                         path_to_output = "data/deconvoluted_epitopes/parsed_deconvoluted_blastout.tsv")
