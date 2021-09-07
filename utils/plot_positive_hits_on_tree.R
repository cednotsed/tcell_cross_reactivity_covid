rm(list=ls())
setwd("../Desktop/git_repos/tcell_cross_reactivity_covid")
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggtree);require(ggplot2);require(seqinr);require(adegenet)
require(ggutils)

# Load data
tree <- read.tree("data/20210803_CoV_genomes_n2572_concat_CGS_aln_trim20_iqtree.delta_rerooted.tree")
epitope.mat <- fread("data/wuhan-hu-1/proteome_blast_heatmap_66_PID.merged.csv")
meta <- fread("data/parsed_host_metadata_FB_LvD_CT.n2572.csv")
mean_df <- fread("results/mean_unexplained_homology_per_accession.csv")
# "EPIISL1098866", "MG772933", 
mean_df <- mean_df %>% 
  mutate(show_point = ifelse(Accession %in% c("MN908947", "EPIISL1098866", 
                                              "MG772933", "MG772934"), 
                             Accession, NA),
         above_t = ifelse(mean > 80, "> 80%", "<= 80%")) %>%
  mutate(show_size = ifelse(!is.na(show_point), 2, 1))
print(all(mean_df$Accession %in% tree$tip.label))

# Match tips to metadata
meta.match <- mean_df[match(tree$tip.label, mean_df$Accession), ]

# Make dataframe for tree annotations
dd <- data.frame(Accession=meta.match$Accession, 
                 show_point=meta.match$show_point,
                 show_size=meta.match$show_size,
                 above_t=meta.match$above_t,
                 mean_hom=meta.match$mean,
                 Genus = meta.match$Genus_clean)
rownames(dd) <- tree$tip.label

# Get tip order to align tips and heatmap
d <- fortify(tree)
d <- subset(d, isTip)
tip_order <- d$label[order(d$y, decreasing = F)]

# Plot positive hits on tree
p <- ggtree(tree)
p1 <- p %<+% dd +
  geom_tippoint(aes(hjust = 0.5,
                    color = show_point,
                    size = show_size), alpha = .8) +
  labs(color = "Viral Genus")  +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  guides(colour = guide_legend(override.aes = list(size = 12)))
# scale_color_manual(values = c("green", "blue","red"), na.value = "grey")

# Plot average HCoV-unexplained PID > 80%
p2 <- p %<+% dd +
  geom_tippoint(aes(hjust = 0.5,
                    color = above_t), alpha = .8) +
  labs(color = "")  +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  guides(colour = guide_legend(override.aes = list(size = 12)))

p2
# ggsave("results/positive_hits_on_tree.png", 
#        width = 10, height = 10, dpi = 300, limitsize = F,
#        plot = p1)

ggsave("results/mean_pid_above_t.png", 
       width = 10, height = 10, dpi = 300, limitsize = F,
       plot = p2)

## INTERACTIVE ##
ggzoom(p)

p_annot <- p + geom_tiplab()

ggzoom(p_annot)


p3 <- p %<+% dd +
  geom_tippoint(aes(hjust = 0.5,
                    color = mean_hom,
                    size = show_size), alpha = .8) +
  labs(color = "Viral Genus")  +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  guides(colour = guide_legend(override.aes = list(size = 12))) +
  scale_color_gradient(low = "blue", high = "red")

p3

ggsave("results/positive_hits_color_by_mean.png", 
       width = 10, height = 10, dpi = 300, limitsize = F,
       plot = p3)
