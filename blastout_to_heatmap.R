# Parse blastout to heatmap matrix 
rm(list=ls())
setwd('~/git_repos/covid_cross_immunity/')
require(tidyverse);require(data.table);require(ape);require(seqinr);require(adegenet)

tree <- read.tree("data/cov_n2531_ml_tree_delta_rerooted.tree")
meta <- fread("data/parsed_host_metadata_FB_LvD_CT.n2531.csv")

blast_to_heatmap <- function(path_to_blastout, path_to_mat, pid_threshold) {
  blastout <- fread(path_to_blastout)
  
  # Match tips to metadata
  meta.match <- meta[match(tree$tip.label, meta$Accession), ]
  
  # Make dataframe for tree annotations
  dd <- data.frame(Accession=meta.match$Accession, 
                   Genus=meta.match$Genus_clean)
  rownames(dd) <- tree$tip.label
  
  # Find unique list of epitopes with hits and construct empty matrix
  epitopes <- unique(blastout$qseqid)
  epitopes <- str_replace(epitopes, "-", "_")
  epitope.mat <- matrix(numeric(), nrow = dim(dd)[1], ncol = length(epitopes)) #empty matrix (columns epitopes, rows accessions)
  colnames(epitope.mat) <- epitopes
  rownames(epitope.mat) <- rownames(dd)
  
  for (i in 1:length(epitopes)) #add epitope counts to empty marix (note choice of pid threshold)
    # for (i in 1:10) #add epitope counts to empty marix (note choice of pid threshold)
  {
    blast.i <- blastout[which(blastout$qseqid == epitopes[i]), ]
    blast.i.perc <- blast.i[which(blast.i$pident > pid_threshold),] #using a pid cut off of 40 but can use whatever
    # print(table(meta[match(blast.i.perc$Accession, meta$Accession), 3])) #provides matches by species assigned in meta
    
    #note not all those which pass the blast filter are in the phylogenetic tree
    #stores only those which match to the phylogeny
    index.match.i <- match(blast.i.perc$Accession, rownames(epitope.mat))
    epitope.mat[na.omit(index.match.i), i] <- blast.i.perc$pident[!is.na(index.match.i)]
  }
  
  epitope.mat <- data.frame(epitope.mat)
  
  # Save data
  fwrite(epitope.mat, path_to_mat, row.names = T)
  
  return(epitope.mat)
}

## Proteome
# proteome_40_mat <- blast_to_heatmap(path_to_blastout = "data/wuhan-hu-1/parsed_proteome_blastout.tsv",
#                                     path_to_mat = "data/wuhan-hu-1/proteome_blast_heatmap_40_PID.merged.csv",
#                                     pid_threshold = 40)
# 
# proteome_66_mat <- blast_to_heatmap(path_to_blastout = "data/wuhan-hu-1/parsed_proteome_blastout.tsv",
#                                     path_to_mat = "data/wuhan-hu-1/proteome_blast_heatmap_66_PID.merged.csv",
#                                     pid_threshold = 66)
# 
# proteome_80_mat <- blast_to_heatmap(path_to_blastout = "data/wuhan-hu-1/parsed_proteome_blastout.tsv",
#                                     path_to_mat = "data/wuhan-hu-1/proteome_blast_heatmap_80_PID.merged.csv",
#                                     pid_threshold = 80)

## 161 deconvoluted epitopes
deconv_40_mat <- blast_to_heatmap(path_to_blastout = "data/deconvoluted_epitopes/parsed_deconvoluted_blastout.tsv",
                                    path_to_mat = "data/deconvoluted_epitopes/deconvoluted_heatmap_40_PID.merged.csv",
                                    pid_threshold = 40)

# deconv_66_mat <- blast_to_heatmap(path_to_blastout = "data/deconvoluted_epitopes/parsed_deconvoluted_blastout.tsv",
#                                     path_to_mat = "data/deconvoluted_epitopes/deconvoluted_heatmap_66_PID.merged.csv",
#                                     pid_threshold = 66)
# 
# deconv_80_mat <- blast_to_heatmap(path_to_blastout = "data/deconvoluted_epitopes/parsed_deconvoluted_blastout.tsv",
#                                     path_to_mat = "data/deconvoluted_epitopes/deconvoluted_heatmap_80_PID.merged.csv",
#                                     pid_threshold = 80)

## Unexposed
# unexposed_40_mat <- blast_to_heatmap(path_to_blastout = "data/unexposed_epitopes/parsed_unexposed_blastout.tsv",
#                                     path_to_mat = "data/unexposed_epitopes/unexposed_blast_heatmap_40_PID.merged.tsv",
#                                     pid_threshold = 40)

# unexposed_66_mat <- blast_to_heatmap(path_to_blastout = "data/wuhan-hu-1/parsed_proteome_blastout.tsv",
#                                     path_to_mat = "data/wuhan-hu-1/proteome_blast_heatmap_66_PID.merged.csv",
#                                     pid_threshold = 66)
# 
# unexposed_80_mat <- blast_to_heatmap(path_to_blastout = "data/wuhan-hu-1/parsed_proteome_blastout.tsv",
#                                     path_to_mat = "data/wuhan-hu-1/proteome_blast_heatmap_80_PID.merged.csv",
#                                     pid_threshold = 80)