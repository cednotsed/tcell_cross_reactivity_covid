## Find out which pools are what ##
rm(list=ls())
setwd('~/git_repos/covid_cross_immunity/')
require(aplot);require(tidyverse);require(reshape2);require(data.table);require(ape);require(ggtree);require(ggplot2);require(seqinr);require(adegenet); require(egg)

# Load data
meta <- fread("data/parsed_host_metadata_FB_LvD.csv")
protein_meta <- fread("data/unexposed_epitopes/unexposed_epitopes.csv")
blastout <- fread("data/unexposed_epitopes/parsed_unexposed_blastout.tsv")


# Get column names of proteins in genome order
pool_order <- c("NSP5", "NSP7", "NSP8", "NSP9", "NSP10", 
                "NSP12-1", "NSP12-2", "NSP12-3", "NSP12-4", "NSP12-5", 
                "NSP13-1", "NSP13-2", "NSP13-3", 
                "NSP14",
                "ALL SPIKE", "ORF3", "M", "ORF7/8", "NP-1", "NP-2")

epitope_order <- c()
for (pool in pool_order) {
  epitope_order <- c(epitope_order, protein_meta$epitope[protein_meta$pool == pool])
}

# Get coordinates for annotating pools on heatmap
len_df <- data.frame(table(protein_meta$pool)[pool_order])
colnames(len_df) <- c("pool", "len")
len_df$cum_len <- 0

# Make len_df cumulative
for (i in seq(nrow(len_df))) {
  len_df[i, "cum_len"] <- sum(len_df[1:i, "len"]) + 0.5
}

plot_df <- blastout %>%
  filter(Accession %in% c("NC_045512", "AY427439")) %>% # Get only SARS2 and SARS1
  select(Accession, qseqid, pident) %>%
  group_by(Accession, qseqid) %>%
  summarise(pident = max(pident, na.rm = T))

# Order heatmap rows and columns
plot_df$qseqid <- factor(plot_df$qseqid, levels = epitope_order)

# Plot Heatmap
hm <- ggplot(plot_df) +
  geom_tile(aes(x = qseqid, y = Accession, fill = pident)) +
  scale_fill_gradient(low = "green", high = "red", na.value = "white") + 
  geom_vline(xintercept = len_df$cum_len, linetype="dotted", color = "black", size = 1.5) +
  theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.key.width=unit(1, "in")) +
  guides(colour = guide_legend(override.aes = list(size = 12))) +
  labs(fill = "Sequence Identity (%)")

hm
ggsave("results/unexposed_epitopes/peptide_to_sars1_2.png", plot = hm, dpi = 600, width = 40, height = 20)

# Check no. of peptides
print(length(unique(blastout$qseqid)))
print(length(unique(plot_df$qseqid)))
print(length(epitope_order))

## Simplify pool homologies #######################################################

# Get total number of epitopes per pool
count_df <- protein_meta %>%
  group_by(pool) %>%
  summarise(expected = n())

# Define pools that are not in both cohorts
not_in_both <- c("NSP5", "NSP8", "NSP9", "NSP10", "NSP14", "ORF7/8")
in_both_df <- count_df %>%
  mutate(Accession = "In both cohorts", 
         perfect = ifelse(!(pool %in% not_in_both), T, F))

# Get list of entries
simple_df <- rbind(count_df %>% mutate(Accession = "NC_045512"),
                   count_df %>% mutate(Accession = "AY427439"))

plot2_df <- plot_df %>% 
  filter(pident == 100) %>%
  rename(epitope = qseqid) %>%
  left_join(protein_meta, by = "epitope") %>%
  group_by(Accession, pool) %>%
  summarise(n_peptides = n()) %>%
  right_join(simple_df, by = c("pool", "Accession")) %>%  # Retain all entries
  mutate(n_peptides = replace_na(n_peptides, 0)) %>%
  mutate(perfect = (n_peptides == expected),
         Accession = recode(Accession, NC_045512 = "SARS-CoV-2", AY427439 = "SARS-CoV-1")) %>%
  select(pool, expected, Accession, perfect)

# plot2_df <- bind_rows(plot2_df, in_both_df)

plot2_df$pool <- factor(plot2_df$pool, levels = pool_order)
in_both_df$pool <- factor(in_both_df$pool, levels = pool_order)
plot2_df$Accession <- factor(plot2_df$Accession, 
                             levels = c("SARS-CoV-1", "SARS-CoV-2", "In both cohorts"))

plt2 <- ggplot(plot2_df, aes(x = pool, y = Accession, fill = perfect)) +
  geom_tile(width = 0.95, height = 1) +
  scale_fill_manual(values = c("grey", "deepskyblue1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_blank()) +
  labs(x = "Peptide Pool", fill = "Perfect homology?")

plt3 <- ggplot(in_both_df, aes(x = pool, y = Accession, fill = perfect)) +
  geom_tile(width = 0.95, height = 1) +
  scale_fill_manual(values = c("grey", "tomato4")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank()) +
  labs(x = "Peptide Pool", fill = "")

plt_final <- egg::ggarrange(plt3, plt2, nrow = 2, heights = c(1/3, 1))

ggsave("results/unexposed_epitopes/pool_homology.png", plot = plt_final, dpi = 300, width = 8, height = 5)

# # epitope blast - read in and make table
# blastout <- fread('data/unexposed_epitopes/blast_unexposed_epitopes.tsv')
# colnames(blastout) <- c('qseqid','sseqid','pident','length','mismatch',
#                         'gapopen', 'qstart', 'qend', 'sstart', 'send',
#                         'evalue', 'bitscore')
# 
# 
# 
# 
# %>%
#   filter(pident == 100) %>%
#   rename(epitope = qseqid) %>%
#   left_join(protein_meta, by = "epitope") %>%
#   group_by(pool) %>%
#   summarise(perfect = n()) %>%
#   right_join(count_df, by = "pool")
# 
# # SARS-CoV-1
# blastout %>%
#   filter(Accession == "AY274119") %>%
#   select(qseqid, pident) %>%
#   group_by(qseqid) %>%
#   summarise(pident = max(pident, na.rm = T)) %>%
#   filter(pident == 100) %>%
#   rename(epitope = qseqid) %>%
#   left_join(protein_meta, by = "epitope") %>%
#   group_by(pool) %>%
#   summarise(perfect = n()) %>%
#   right_join(count_df, by = "pool")