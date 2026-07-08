cd /home/matiasmedina/salmon_multiomics/results/02_fastp/

for i in {01..23}; do
    echo "ALINEANDO MUESTRA: $i"
    
    STAR --runThreadN 8 \
         --genomeDir /home/matiasmedina/salmon_multiomics/data/reference/star_index/ \
         --readFilesIn Muestra_${i}_R1_clean.fastq.gz Muestra_${i}_R2_clean.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix /home/matiasmedina/salmon_multiomics/results/05_aligned_bams/Muestra_${i}_ \
         --outSAMtype BAM SortedByCoordinate \
         --twopassMode Basic \
         --quantMode GeneCounts
done