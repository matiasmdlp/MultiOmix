# 01: Análisis Diferencial de Expresión (Transcriptómica)

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(ggvenn)
library(patchwork)

# 1. Carga y Sincronización de Datos
metadata <- read.csv2("data/metadata.csv", row.names = 1, check.names = FALSE)

# Convertir columnas categóricas a factores
cols_to_factor <- c("Regimen", "Empresa", "Orden", "Region", "Centro", 
                    "Temporada", "Mediciones", "Organismo", "Estado", 
                    "Primera_Medicion_Innovex", "Ultima_Medicion_Innovex", 
                    "Primera_Fecha_Seleccionada", "Ultima_Fecha_Seleccionada")
metadata[cols_to_factor] <- lapply(metadata[cols_to_factor], as.factor)

raw_data <- read.table("data/salmon_counts.txt", header = TRUE, skip = 1, row.names = 1, check.names = FALSE)
colnames(raw_data) <- basename(colnames(raw_data))
colnames(raw_data) <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam", "", colnames(raw_data))

# Sincronización
muestras_comunes <- intersect(rownames(metadata), colnames(raw_data))
counts_final <- as.matrix(raw_data[, muestras_comunes])
mode(counts_final) <- "integer"
metadata <- metadata[muestras_comunes, ]

# 2. DESeq2 
dds <- DESeqDataSetFromMatrix(countData = counts_final, colData = metadata, design = ~ Regimen + Temporada)
dds <- DESeq(dds)
save(dds, metadata, file = "results/checkpoint_final.RData")

# 3. PCA Global
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Temporada", "Regimen"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Temporada, shape=Regimen)) +
  geom_point(size=4, alpha=0.8) +
  coord_fixed(ratio = 1, xlim = c(-65, 65), ylim = c(-65, 65)) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") + 
  labs(title = "Análisis de Componentes Principales (PCA)",
       x = paste0("PC1: ", percentVar[1], "% varianza"),
       y = paste0("PC2: ", percentVar[2], "% varianza")) +
  theme(aspect.ratio = 1, legend.position = "right")
ggsave("results/plots/PCA_Estandarizado_Final.png", width = 8, height = 8, dpi = 300)

# 4. Función de Volcano Plots
generar_reporte_volcano <- function(res_obj, nombre_grafico) {
  message(paste("Generando gráfico para:", nombre_grafico))
  write.csv(as.data.frame(res_obj), file = paste0("results/tables/Tabla_", nombre_grafico, ".csv"))
  
  p <- EnhancedVolcano(res_obj,
    lab = rownames(res_obj), x = 'log2FoldChange', y = 'padj',
    title = nombre_grafico, pCutoff = 0.05, FCcutoff = 2,
    xlim = c(-12, 12), ylim = c(0, 15),
    pointSize = 1.5, labSize = 3.0, drawConnectors = TRUE)
  
  ggsave(filename = paste0("results/plots/Volcano_", nombre_grafico, ".png"), plot = p, width = 10, height = 8)
}

# Ejecución para Temporadas
niveles <- levels(metadata$Temporada)
combos <- combn(niveles, 2)
for(i in 1:ncol(combos)){
  nivel_A <- combos[1, i]
  nivel_B <- combos[2, i]
  id_nombre <- paste0("Temporada_", nivel_A, "_vs_", nivel_B)
  res_actual <- results(dds, contrast = c("Temporada", nivel_A, nivel_B))
  generar_reporte_volcano(res_actual, id_nombre)
}