## Script for plotting ML tree and blast output heatmap ##
rm(list=ls())
setwd("~/git_repos/tcell_cross_reactivity_covid")
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggtree);require(ggplot2);require(seqinr);require(adegenet)

# Load data
tree <- read.tree("data/cov_n2531_ml_tree_delta_rerooted.tree")
epitope.mat <- fread("data/deconvoluted_epitopes/mateus_nina_nelde_blast_heatmap_40_PID.merged.csv")
meta <- fread("data/parsed_host_metadata_FB_LvD_CT.n2531.csv")

# Match tips to metadata
meta.match <- meta[match(tree$tip.label, meta$Accession), ]

# Make dataframe for tree annotations
dd <- data.frame(Accession=meta.match$Accession, 
                 Genus=meta.match$Genus_clean)
rownames(dd) <- tree$tip.label

# Read data and convert to long format
plot_df <- epitope.mat %>% 
  gather(key = Accession, value = value, -V1) %>%
  rename(Accession = V1, Query = Accession) %>%
  as_tibble()

plot_df <- plot_df %>%
  separate(Query, into = c(NA, "protein", NA), sep = "_", remove = F)

# Get coordinates for annotating pools on heatmap
protein_order <- c("NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP10", 
                   "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", 
                   "S", "ORF3", "E", "M", "ORF6", "ORF7", "ORF8", "N")

epitope_order <- c()
len <- c()

for (p in protein_order) {
  epitopes <- unique(plot_df$Query[plot_df$protein == p])
  epitope_order <- c(epitope_order, epitopes)
  len <- c(len, length(epitopes))
}

len_df <- data.frame(protein = protein_order, len = len)
len_df$cum_len <- 0

# Make len_df cumulative
for (i in seq(nrow(len_df))) {
  len_df[i, "cum_len"] <- sum(len_df[1:i, "len"]) + 0.5
}

# Get tip order to align tips and heatmap
d <- fortify(tree)
d <- subset(d, isTip)
tip_order <- d$label[order(d$y, decreasing = F)]

# Order heatmap rows, columns, grids
plot_df$Accession <- factor(plot_df$Accession, levels = tip_order)
plot_df$protein <- factor(plot_df$protein, levels = protein_order)
plot_df$Query <- factor(plot_df$Query, levels = epitope_order)

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
  geom_vline(xintercept = len_df$cum_len, linetype="dashed", color = "black", size = 1.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5),
        axis.title = element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.key.width=unit(1, "in")) +
  guides(colour = guide_legend(override.aes = list(size = 12))) +
  labs(fill = "Sequence Identity (%)")

hm_with_annot <- ggplot(plot_df) +
  geom_tile(aes(x = Query, y = Accession, fill = value)) +
  scale_fill_gradient(low = "green", high = "red", na.value = "white") + 
  geom_vline(xintercept = len_df$cum_len, linetype="dashed", color = "black", size = 1.5) +
  theme(axis.text.y = element_blank(),
        # axis.ticks = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5),
        axis.title = element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.key.width=unit(1, "in")) +
  guides(colour = guide_legend(override.aes = list(size = 12))) +
  labs(fill = "Sequence Identity (%)")


ggsave("results/deconvoluted_epitopes/hm_with_text.merged.png", plot = hm_with_annot, dpi = 200, width = 15, height = 15)

# Combine plots
png('results/deconvoluted_epitopes/mateus_nina_nelde_sars_2_genome_heatmap_raxml_delta_rooted_40_PID.merged.png', width = 40, height = 15, units = 'in', res = 300)
par(mar = c(10.1, 4.1, 4.1, 2.1), xpd = T)
hm %>% insert_left(p)
dev.off()
