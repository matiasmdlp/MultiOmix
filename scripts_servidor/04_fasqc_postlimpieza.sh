fastqc -o /home/matiasmedina/salmon_multiomics/results/04_fastqc_clean/ \
       /home/matiasmedina/salmon_multiomics/results/02_fastp/*.fastq.gz

cd /home/matiasmedina/salmon_multiomics/
multiqc results/ -o results/03_multiqc/ -f