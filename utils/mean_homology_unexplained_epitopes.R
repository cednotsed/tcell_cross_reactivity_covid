setwd("../Desktop/git_repos/tcell_cross_reactivity_covid")
require(ggplot2); require(tidyverse); require(data.table)
epitope_df <- fread("data/deconvoluted_epitopes/parsed_deconvoluted_blastout.tsv")
unexplained_epitopes <- fread("data/deconvoluted_epitopes/list_of_unexplained_epitopes.txt", header = F)$V1

viruses <- c("MN996532", "EPIISL410721", "EPIISL412977",
             "EPIISL1098866", "MG772933", "MG772934")

pos_hits <- c("EPIISL1098866", "MG772933", "MG772934")

hit_annot <- epitope_df %>%
  filter(qseqid %in% unexplained_epitopes) %>%
  group_by(Accession) %>%
  summarise(mean = mean(pident)) %>%
  filter(Accession %in% viruses) %>%
  mutate(Accession = case_when(
    grepl("MN996532", Accession) ~ "RaTG13",
    grepl("EPIISL410721", Accession) ~ "Pangolin-CoV",
    grepl("EPIISL412977", Accession) ~ "RmYN02",
    grepl("EPIISL1098866", Accession) ~ "PrC31",
    grepl("MG772933", Accession) ~ "Bat-SL-CoVZC45",
    grepl("MG772934", Accession) ~ "Bat-SL-CoVZXC21"))

hit_line <- hit_annot$mean

# Mean unexplained homology per accession
mean_df <- epitope_df %>% 
  group_by(Accession) %>%
  summarise(mean = mean(pident))

mean_meta <- mean_df %>% 
  left_join(fread("data/parsed_host_metadata_FB_LvD_CT.n2572.csv"))
# fwrite(mean_meta, "results/mean_unexplained_homology_per_accession.csv")


h <- mean_df %>%
  ggplot(aes(x = mean)) +
  geom_histogram(bins = 500) +
  labs(x = "Mean homology of unexplained epitopes", y = "Frequency") +
  geom_vline(xintercept = hit_line, color = "red", lty = "dotted") +
  geom_text(data = hit_annot, 
            aes(x = mean, y = 100, label = Accession), 
            color = "red",
            size = 2,
            position = position_jitter(height = 20))
h

ggsave("results/mean_unexplained_epitope_homology.png", dpi = 300, plot = h, width = 8, height = 5)

hookup_table <- fread("data/SARS-CoV-2_hookup_table.csv") %>%
  distinct(protein_name, codon_number, .keep_all = T)

epitope_merged <- tibble(unexplained_epitopes) %>%
  separate("unexplained_epitopes", into = c("study", "protein_name", "codon_number"), sep = "_") %>%
  mutate(codon_number = as.numeric(codon_number)) %>%
  left_join(hookup_table, by = c("protein_name", "codon_number"))

# Calculate average homology of epitope hits
epitope_df %>%
  filter(qseqid %in% unexplained_epitopes, Accession %in% viruses) %>%
  separate("qseqid", 
           into = c("study", "protein_name", "codon_number"), 
           sep = "_", 
           remove = F) %>%
  filter(protein_name %in% c("S", "N")) %>%
  mutate(Accession = case_when(
    grepl("MN996532", Accession) ~ "RaTG13",
    grepl("EPIISL410721", Accession) ~ "Pangolin-CoV",
    grepl("EPIISL412977", Accession) ~ "RmYN02",
    grepl("EPIISL1098866", Accession) ~ "PrC31",
    grepl("MG772933", Accession) ~ "Bat-SL-CoVZC45",
    grepl("MG772934", Accession) ~ "Bat-SL-CoVZXC21")) %>%
  group_by(Accession) %>%
  summarise(mean = mean(pident))


# Mean 
epitope_df %>%
  filter(qseqid %in% unexplained_epitopes, Accession %in% pos_hits) %>%
  separate(qseqid, into = c("study", "protein", NA), sep = "_", remove = F) %>%
  group_by(Accession) %>%
  mutate(protein = ifelse(grepl("NSP", protein), "ORF1ab", protein),
         above_t = ifelse(pident > 80, "PID > 80%", "PID <= 80%")) %>%
  count(protein, above_t) %>%
  ggplot(aes(x = Accession, y = n, fill = protein)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(above_t)) +
  labs(y = "No. of unexplained epitopes", fill = "Protein") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Mean pident per epitope
order_df <- epitope_df %>%
  separate(qseqid, into = c("study", "protein", NA), sep = "_", remove = F)
qseqid_order <- unique(order_df$qseqid[order(order_df$protein)])

epitope_df %>%
  filter(qseqid %in% unexplained_epitopes) %>%
  mutate(qseqid = factor(qseqid, levels = qseqid_order)) %>%
  group_by(qseqid) %>%
  summarise(mean = mean(pident), min = min(pident), max = max(pident)) %>%
  ggplot(aes(x = qseqid, y = mean)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot counts of unexplained epitopes above and below 80% homology
unexplained_plt <- epitope_df %>%
  filter(qseqid %in% unexplained_epitopes) %>%
  left_join(fread("data/parsed_host_metadata_FB_LvD_CT.n2572.csv")) %>%
  mutate(above_t = ifelse(pident > 80, "PID > 80%", "PID <= 80%"),
         qseqid = factor(qseqid, levels = qseqid_order)) %>%
  group_by(qseqid, above_t, Genus_clean) %>%
  count(above_t) %>%
  ggplot(aes(x = qseqid, y = n, fill = Genus_clean)) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(above_t)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Epitope", y = "No. of Accessions", title = "Unexplained epitopes") +
  ylim(0, 2500)

explained_plt <- epitope_df %>%
  filter(!(qseqid %in% unexplained_epitopes)) %>%
  left_join(fread("data/parsed_host_metadata_FB_LvD_CT.n2572.csv")) %>%
  mutate(above_t = ifelse(pident > 80, "PID > 80%", "PID <= 80%"),
         qseqid = factor(qseqid, levels = qseqid_order)) %>%
  group_by(qseqid, above_t, Genus_clean) %>%
  count(above_t) %>%
  ggplot(aes(x = qseqid, y = n, fill = Genus_clean)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(above_t)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Epitope", y = "No. of Accessions", title = "Explained epitopes") +
  ylim(0, 2500)

plt <- ggpubr::ggarrange(unexplained_plt, explained_plt)

ggsave("results/num_accessions_per_epitope.png", dpi = 300, width = 30, height = 10)
