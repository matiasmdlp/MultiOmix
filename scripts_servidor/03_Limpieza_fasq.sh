cd /home/matiasmedina/salmon_multiomics/data/raw_renamed/

for i in {01..23}; do
    echo "PROCESANDO FILTROS MUESTRA: $i"
    
    fastp \
        -i Muestra_${i}_R1.fastq.gz \
        -I Muestra_${i}_R2.fastq.gz \
        -o /home/matiasmedina/salmon_multiomics/results/02_fastp/Muestra_${i}_R1_clean.fastq.gz \
        -O /home/matiasmedina/salmon_multiomics/results/02_fastp/Muestra_${i}_R2_clean.fastq.gz \
        --trim_poly_g \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 50 \
        -h /home/matiasmedina/salmon_multiomics/results/02_fastp/Muestra_${i}.html \
        -j /home/matiasmedina/salmon_multiomics/results/02_fastp/Muestra_${i}.json \
        --thread 4
done