## Script for plotting ML tree and blast output heatmap ##
rm(list=ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggtree);require(ggplot2);require(seqinr);require(adegenet)

# Load data
tree <- read.tree("data/cov_n2531_ml_tree_delta_rerooted.tree")
epitope.mat <- fread("data/wuhan-hu-1/proteome_blast_heatmap_80_PID.merged.csv")
meta <- fread("data/parsed_host_metadata_FB_LvD_CT.n2531.csv")

print(all(meta$Accession %in% tree$tip.label))

# Match tips to metadata
meta.match <- meta[match(tree$tip.label, meta$Accession), ]

# Make dataframe for tree annotations
dd <- data.frame(Accession=meta.match$Accession, 
                 Genus=meta.match$Genus_clean)
rownames(dd) <- tree$tip.label

# Get column names of proteins in genome order
protein_vec <- fread("data/wuhan-hu-1/wuhan-hu-1_genome_annotations.csv")$Name

colname_vec <- c()
len_df <- data.frame()

for (prot_name in protein_vec) {
  to_append <- colnames(epitope.mat)[grepl(paste0(prot_name, "_"), colnames(epitope.mat))]
  colname_vec <- c(colname_vec, to_append)
  len_df <- rbind(len_df, data.frame(protein = c(prot_name), len = c(length(colname_vec))))
}

# Read data and convert to long format
plot_df <- epitope.mat %>% 
  gather(key = Accession, value = value, -V1) %>%
  rename(Accession = V1, Query = Accession) %>%
  as.data.frame()

# Get tip order to align tips and heatmap
d <- fortify(tree)
d <- subset(d, isTip)
tip_order <- d$label[order(d$y, decreasing = F)]

# Order heatmap rows and columns
plot_df$Accession <- factor(plot_df$Accession, levels = tip_order)
plot_df$Query <- factor(plot_df$Query, levels = colname_vec)

# Plot tree
p <- ggtree(tree)
p <- p %<+% dd +
  geom_tippoint(aes(hjust = 0.5,
                    color = Genus), alpha = .6) +
  labs(color = "Viral Genus")  +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  guides(colour = guide_legend(override.aes = list(size = 12)))

# Plot Heatmap
hm <- ggplot(plot_df) +
  geom_tile(aes(x = Query, y = Accession, fill = value)) +
  scale_fill_gradient(low = "green", high = "red", na.value = "white") + 
  geom_vline(xintercept = len_df$len + 0.5, linetype="dotted", color = "black", size = 1.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.key.width=unit(1, "in")) +
  guides(colour = guide_legend(override.aes = list(size = 12))) +
  labs(fill = "Sequence Identity (%)")

# Combine plots
png("results/proteome_heatmap_raxml_rooted_80_PID.merged.png", width = 40, height = 15, units = 'in', res = 300)
par(mar = c(10.1, 4.1, 4.1, 2.1), xpd = T)
hm %>% insert_left(p)
dev.off()

########################## TEST ###################################################################
# test_tree <- rtree(15)
# 
# test_meta <- data.frame(label = test_tree$tip.label)
# test_meta$Genus <- c(rep("Mus", 5), rep("Homo", 5), rep("Xenopus", 5))
# 
# p_test <- ggtree(test_tree) %<+% test_meta +
#   geom_tippoint(aes(color = Genus), size = 4) +
#   geom_tiplab(hjust=-0.5) +
#   theme(legend.title = element_text(size = 30),
#         legend.text = element_text(size = 30)) +
#   guides(colour = guide_legend(override.aes = list(size=12)))
# 
# hm <- plot_df[1:length(test_tree$tip.label), ]  
# hm$Accession <- test_tree$tip.label
# d <- fortify(test_tree)
# d <- subset(d, isTip)
# tip_order <- d$label[order(d$y, decreasing=F)]
# 
# hm$Accession <- factor(hm$Accession, levels = tip_order)
# 
# test_hm <- ggplot(hm) +
#   geom_tile(aes(x = Query, y = Accession, fill = value)) +
#   theme(legend.title = element_text(size = 30),
#         legend.text = element_text(size = 25),
#         legend.key.width=unit(1, "in")) +
#   guides(colour = guide_legend(override.aes = list(size=30))) +
#   labs(fill = "Sequence Identity (%)")
# 
# test_hm
# pdf('results/test2.pdf',40,15)
# test_hm %>% insert_left(p_test)
# dev.off()
