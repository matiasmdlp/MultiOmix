# MultiOmix
# Pipeline Multiómico para Salmo salar 

Este repositorio contiene el código fuente y el flujo de trabajo analítico desarrollado como parte de la memoria de título para optar al título de Ingeniero Civil Informático. El proyecto se centra en el procesamiento bioinformático y la integración de datos transcriptómicos y proteómicos de *Salmo salar*.

El flujo de trabajo está dividido en dos grandes etapas: el preprocesamiento de lecturas crudas en un entorno de servidor (Bash) y el análisis estadístico multiómico local (R).

## 1. Flujo de Trabajo en Servidor (Preprocesamiento Transcriptómico)
Ubicados en la carpeta `scripts_servidor/`, estos scripts están diseñados para ejecutarse de manera secuencial en un servidor o cluster (ej. servidor `ada`):

* **01_crear_carpetas.sh**: Construye la estructura de directorios necesaria para almacenar los resultados intermedios y el índice de STAR.
* **02_fastqc_raw.sh**: Ejecuta el control de calidad inicial sobre las lecturas `.fastq.gz` crudas.
* **03_Limpieza_fasq.sh**: Realiza el recorte de adaptadores y filtrado por calidad (trimming) utilizando `fastp`.
* **04_fasqc_postlimpieza.sh**: Verifica la calidad de las lecturas limpias y consolida los reportes utilizando `MultiQC`.
* **05_indice_star.sh**: Genera el índice del genoma de referencia de *Salmo salar* (GCF_905237065.1_Ssal_v3.1) optimizando el uso de RAM.
* **06_alineamiento_star.sh**: Alinea las lecturas limpias contra el genoma de referencia, generando archivos BAM ordenados por coordenadas.
* **07_featurecounts.sh**: Cuantifica la expresión génica a nivel de exón, generando la matriz de conteos final (`matriz_counts_final.txt`) para el análisis en R.

*(Nota: Para descargar la matriz de conteos al equipo local, se puede utilizar el comando `scp` incluido en los comentarios del script 07).*

## 2. Flujo de Trabajo Local (Análisis Diferencial e Integración)
Ubicados en la carpeta `scripts_R/`, estos scripts toman la matriz de conteos y los datos proteómicos para el modelado biológico:

* **01_transcriptomica_deseq2.R**: Implementación de `DESeq2` para la normalización de conteos de ARN, análisis de componentes principales (PCA) y evaluación de expresión diferencial estacional mediante Volcano Plots.
* **02_integracion_mixomics.R**: Manejo de la matriz proteómica (imputación de valores perdidos por límite de detección) e integración hologenómica utilizando el framework DIABLO (`mixOmics`) para encontrar firmas altamente correlacionadas.

## Estructura de Datos (No versionada)
Debido a restricciones de tamaño, los siguientes directorios no se suben al repositorio (manejados vía `.gitignore`):
* `data/`: Contiene los archivos crudos (FASTQ), el genoma de referencia (FASTA/GTF) y la metadata.
* `results/`: Almacena todos los archivos BAM, reportes HTML pesados y salidas intermedias generadas por los scripts del servidor.
* `results_R`: [Drive](https://drive.google.com/drive/folders/1xrE7PHYMDAAVpgiaCY7WSMcxrArmTxX0?usp=sharing)

## Requisitos y Dependencias
**Entorno de Servidor (Linux/Bash):**
* FastQC
* fastp
* MultiQC
* STAR
* Subread (featureCounts)

**Entorno Local (R >= 4.0.0):**
* **Bioconductor:** `DESeq2`, `mixOmics`, `EnhancedVolcano`, `SummarizedExperiment`.
* **CRAN:** `ggplot2`, `pheatmap`, `ggvenn`, `readxl`.
