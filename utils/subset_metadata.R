rm(list = ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(tidyverse); require(data.table); require(ggplot2); require(seqinr); require("ape")
meta <- fread("data/parsed_host_metadata_FB_LvD.csv")
tree <- read.tree("data/cov_n2531_ml_tree_delta_rerooted.tree")

meta_filt <- meta %>%
  filter(Accession %in% tree$tip.label) %>%
  select(-Release_Date, -Family, -Sequence_Type, 
         -Segment, -Nuc_Completeness, -y,
         -BioSample, -Genotype)


fwrite(meta_filt, "data/parsed_host_metadata_FB_LvD_CT.n2531.csv")
