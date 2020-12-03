awk 'BEGIN { FS = "," }; NR>1{if ($2 != "") {print $2}}' wuhan-hu-1_genome_annotations.csv > wuhan-hu-1_protein_accessions.txt
