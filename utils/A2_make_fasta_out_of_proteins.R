# Parse individual Wuhan-Hu-1 proteins into 15-mer peptides
rm(list = ls())
setwd("../data/wuhan-hu-1")
require(data.table)
require(tidyverse)
require(Biostrings)

kmer_length = 15

fasta_vec <- c()
for (file in list.files("individual_proteins/")) {
  print(file)
  fasta <- read.fasta(paste0("individual_proteins/", file), as.string = F, set.attributes = F)[[1]]
  sequence <- toupper(fasta)
  
  if (length(sequence) >= 15) {
    for(i in seq(length(sequence) - (kmer_length - 1))) {
      kmer_name <- paste0(">", str_remove(file, ".fasta"), "_", i, "_", i + 14)
      kmer <- paste(sequence[i: (i + 14)], collapse="")
      fasta_vec <- c(fasta_vec, kmer_name, kmer)
    }
  }
}

writeLines(fasta_vec, "wuhan-hu-1_epitopes.fasta")
