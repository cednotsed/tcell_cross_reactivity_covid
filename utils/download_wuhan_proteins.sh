base_dir=../data/wuhan-hu-1

for line in $(awk 'BEGIN { FS = "," }; NR>1 {if ($2 != "") {print}}' ${base_dir}/wuhan-hu-1_genome_annotations.csv)
do
    protein=$(echo $line| awk 'BEGIN { FS = "," }; {print $1}')
    acc=$(echo $line| awk 'BEGIN { FS = "," }; {print $2}')
    echo $acc
    esearch -db protein -query ${acc}| efetch -format fasta > ${base_dir}/${protein}.fasta
done
