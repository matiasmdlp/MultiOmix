cd /home/matiasmedina/salmon_multiomics/results/06_counts

featureCounts \
  -T 8 \
  -p \
  -s 0 \
  -M \
  --fraction \
  -t exon \
  -g gene \
  -a /home/matiasmedina/salmon_multiomics/data/reference/GCF_905237065.1_Ssal_v3.1_genomic.gtf \
  -o matriz_counts_final.txt \
  /home/matiasmedina/salmon_multiomics/results/05_aligned_bams/*.bam


# Si necesitas descar los resultados, usa el comando en la terminal local de tu equipo (recuerda cambiar el usaurio y carpetas)
# scp matiasmedina@ada:/home/matiasmedina/salmon_multiomics/results/06_counts/matriz_counts_final.txt ~/Downloads/