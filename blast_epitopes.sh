## Make db ##
#makeblastdb -in data/entrez_proteins/protein_seqs_n2804_120620.fasta -parse_seqids -title "CoV Proteome" -out data/entrez_proteins/cov_protein_2804_blast_db/protein_2804_blast_db -dbtype prot

#makeblastdb -in data/entrez_proteins/protein_seqs_n2531.fasta -parse_seqids -title "CoV Proteome" -out data/entrez_proteins/cov_protein_2531_blast_db/protein_2531_blast_db -dbtype prot

#makeblastdb -in data/ncbi_virus_proteins/coronaviridae_ncbi_virus_proteins_231120.fasta -parse_seqids -title "NCBI virus" -out data/ncbi_virus_proteins/NCBI_virus_CoV_protein_blast_db/NCBI_virus_CoV_protein_blast_db -dbtype prot

#makeblastdb -in data/prokka_proteins/prokka_proteins_2531.fasta -parse_seqids -title "Prokka" -out data/prokka_proteins/prokka_2531_CoV_protein_blast_db/prokka_2531_CoV_protein_blast_db -dbtype prot

#makeblastdb -in data/SARS1_SARS2.fasta -parse_seqids -title "SARS1_2 Proteome" -out data/SARS1_2_blast_db/SARS1_2_blast_db -dbtype prot

## Blastp SARS1 and SARS2 ##
#db=data/SARS1_2_blast_db/SARS1_2_blast_db
#out=data/unexposed_epitopes/blast_sars1_sars2.prokka_db.tsv
#query=data/unexposed_epitopes/unexposed_epitopes.fasta

## Blastp batch entrez annotations ##
#db=data/entrez_proteins/cov_protein_2531_blast_db/protein_2531_blast_db

#out=data/wuhan-hu-1/blast_wuhan_epitopes.batch_entrez.tsv
#query=data/wuhan-hu-1/wuhan-hu-1_epitopes.fasta

#out=data/deconvoluted_epitopes/blast_deconvoluted_epitopes.batch_entrez.tsv
#query=data/deconvoluted_epitopes/deconvoluted_epitopes.fasta

## Blastp prokka annotations ##
db=data/prokka_proteins/prokka_2531_CoV_protein_blast_db/prokka_2531_CoV_protein_blast_db

#out=data/wuhan-hu-1/blast_wuhan_epitopes.prokka_db.tsv
#query=data/wuhan-hu-1/wuhan-hu-1_epitopes.fasta

out=data/deconvoluted_epitopes/blast_deconvoluted_epitopes.prokka_db.tsv
query=data/deconvoluted_epitopes/deconvoluted_epitopes.fasta

#~/ncbi-blast-2.11.0+/bin/blastp -db ${db} \
#	-task blastp \
#    -query ${query} \
#    -out ${out} \
#    -outfmt 6 \
#    -evalue 200000 \
#	-word_size 2 \
#	-threshold 11 \
#	-window_size 40 \
#	-gapopen 9 \
#	-gapextend 1 \
#	-matrix PAM30 \
#	-comp_based_stats 0 \
#    -qcov_hsp_perc 99 \
#    -num_alignments 1000000000 \
#	-num_threads 12

~/ncbi-blast-2.11.0+/bin/blastp -db ${db} \
	-task blastp-short \
   -query ${query} \
   -out ${out} \
    -outfmt 6 \
    -evalue 9999999999 \
    -qcov_hsp_perc 99 \
    -num_alignments 1000000000 \
    -num_threads 12

