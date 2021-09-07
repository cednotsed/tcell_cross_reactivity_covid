## Make db ##
#makeblastdb -in data/entrez_proteins/protein_seqs_n2572.fasta -parse_seqids -title "CoV Proteome" -out data/entrez_proteins/cov_protein_2572_blast_db/protein_2572_blast_db -dbtype prot

#makeblastdb -in data/prokka_proteins/prokka_proteins_2572.fasta -parse_seqids -title "Prokka" -out data/prokka_proteins/prokka_2572_CoV_protein_blast_db/prokka_2572_CoV_protein_blast_db -dbtype prot

#makeblastdb -in data/SARS1_SARS2.fasta -parse_seqids -title "SARS1_2 Proteome" -out data/SARS1_2_blast_db/SARS1_2_blast_db -dbtype prot

## Blastp batch entrez annotations ##
#db=data/entrez_proteins/cov_protein_2572_blast_db/protein_2572_blast_db

#out=data/wuhan-hu-1/blast_wuhan_epitopes.batch_entrez.tsv
#query=data/wuhan-hu-1/wuhan-hu-1_epitopes.fasta

#out=data/deconvoluted_epitopes/blast_deconvoluted_epitopes.batch_entrez.tsv
#query=data/deconvoluted_epitopes/deconvoluted_epitopes.fasta

## Blastp prokka annotations ##
db=data/prokka_proteins/prokka_2572_CoV_protein_blast_db/prokka_2572_CoV_protein_blast_db

#out=data/wuhan-hu-1/blast_wuhan_epitopes.prokka_db.tsv
#query=data/wuhan-hu-1/wuhan-hu-1_epitopes.fasta

out=data/deconvoluted_epitopes/blast_deconvoluted_epitopes.prokka_db.tsv
query=data/deconvoluted_epitopes/deconvoluted_epitopes.fasta

echo ${query}
echo ${out}
echo ${db}

~/ncbi-blast-2.11.0+/bin/blastp -db ${db} \
	-task blastp-short \
   -query ${query} \
   -out ${out} \
    -outfmt 6 \
    -evalue 9999999999 \
    -qcov_hsp_perc 99 \
    -num_alignments 1000000000 \
    -num_threads 12

