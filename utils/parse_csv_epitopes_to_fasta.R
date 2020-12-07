rm(list = ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(tidyverse); require(data.table); require(ggplot2); require(seqinr)

deconv_df <- fread("data/deconvoluted_epitopes/deconvoluted_epitopes.csv")

deconv_df <- deconv_df %>% 
  unite("query", study, protein, start, sep = "_")

write.fasta(as.list(deconv_df$sequence), deconv_df$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.fasta")

# Split into subsets for blasting online
deconv_df1 <- deconv_df[1:20, ]
deconv_df2 <- deconv_df[21:40, ]
deconv_df3 <- deconv_df[41:60, ]
deconv_df4 <- deconv_df[61:80, ]
deconv_df5 <- deconv_df[81:100, ]
deconv_df6 <- deconv_df[101:120, ]
deconv_df7 <- deconv_df[121:140, ]
deconv_df8 <- deconv_df[141:160, ]
deconv_df9 <- deconv_df[161:175, ]

write.fasta(as.list(deconv_df1$sequence), deconv_df1$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part1.fasta")
write.fasta(as.list(deconv_df2$sequence), deconv_df2$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part2.fasta")
write.fasta(as.list(deconv_df3$sequence), deconv_df3$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part3.fasta")
write.fasta(as.list(deconv_df4$sequence), deconv_df4$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part4.fasta")
write.fasta(as.list(deconv_df5$sequence), deconv_df5$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part5.fasta")
write.fasta(as.list(deconv_df6$sequence), deconv_df6$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part6.fasta")
write.fasta(as.list(deconv_df7$sequence), deconv_df7$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part7.fasta")
write.fasta(as.list(deconv_df8$sequence), deconv_df8$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part8.fasta")
write.fasta(as.list(deconv_df9$sequence), deconv_df9$query, "data/deconvoluted_epitopes/deconvoluted_epitopes.part9.fasta")