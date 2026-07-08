cd /home/matiasmedina/salmon_multiomics/data/reference/

STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir ./star_index \
     --genomeFastaFiles GCF_905237065.1_Ssal_v3.1_genomic.fna \
     --sjdbGTFfile GCF_905237065.1_Ssal_v3.1_genomic.gtf \
     --sjdbOverhang 150 \
     --genomeSAsparseD 2 \
     --limitGenomeGenerateRAM 28000000000