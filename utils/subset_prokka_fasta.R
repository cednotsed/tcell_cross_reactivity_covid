rm(list = ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(tidyverse); require(data.table); require(ggplot2); require(seqinr); require("ape")
meta <- fread("data/parsed_host_metadata_FB_LvD.csv")
tree <- read.tree("data/cov_n2531_ml_tree_delta_rerooted.tree")
protein <- read.fasta("data/prokka_proteins/prokka_proteins_2544.fasta")
prot_names <- names(protein)

NC <- prot_names[grepl("NC_", prot_names)]
rest <- prot_names[!grepl("NC_", prot_names)]
NC_acc <- colsplit(NC, "_", c("desired1", "desired2", NA))[, 1:2]
NC_acc <- paste0(NC_acc$desired1, "_0", NC_acc$desired2)
rest_acc <- colsplit(rest, "_", c("desired", NA))$desired
acc <- unique(c(NC_acc, rest_acc))
to_remove <- acc[!(acc %in% tree$tip.label)]
regex <- paste0(to_remove, collapse = "|")
prot_names_filt <- prot_names[!grepl(regex, prot_names)]

prot_names_filt
prot_filt <- protein[prot_names %in% prot_names_filt]
write.fasta(prot_filt, prot_names_filt, "data/prokka_proteins/prokka_proteins_2531.fasta")
