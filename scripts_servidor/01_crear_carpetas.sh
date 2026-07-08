cd /home/matiasmedina/salmon_multiomics/

# Crear la estructura de directorios para los resultados
mkdir -p results/01_fastqc_raw \
         results/02_fastp \
         results/03_multiqc \
         results/04_fastqc_clean \
         results/05_aligned_bams \
         results/06_counts

# Crear la carpeta para el índice de STAR 
mkdir -p data/reference/star_index