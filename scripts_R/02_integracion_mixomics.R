
# 02: Preprocesamiento Proteómico e Integración (MixOmics)

library(readxl)
library(mixOmics)
library(SummarizedExperiment)

# 1. Carga y Preprocesamiento 
prot_data <- read_excel("data/matrix_proteomico.xlsx")
prot_raw <- as.data.frame(prot_data)
rownames(prot_raw) <- prot_raw[, 1]

columnas_muestras <- grep("Muestra", colnames(prot_raw))
prot_matrix <- prot_raw[, columnas_muestras]
prot_matrix <- prot_matrix[, rownames(metadata)] 

# Filtrado e Imputación
porcentaje_NAs <- apply(prot_matrix, 1, function(x) sum(is.na(x)) / length(x))
prot_filtrada <- prot_matrix[porcentaje_NAs < 0.30, ]
prot_log2 <- log2(prot_filtrada)

X_prot_final <- apply(prot_log2, 1, function(x) {
  minimo_prot <- min(x, na.rm = TRUE)
  x[is.na(x)] <- minimo_prot
  return(x)
})

# 2. Preparación Bloques Ómicos
load("results/checkpoint_final.RData") 
vsd <- vst(dds, blind = FALSE)
varianzas <- apply(assay(vsd), 1, var)
top_genes <- names(sort(varianzas, decreasing = TRUE)[1:3000])
X_rna <- t(assay(vsd)[top_genes, ])

X_prot_final <- X_prot_final[rownames(X_rna), ]

X <- list(mRNA = X_rna, Proteinas = X_prot_final)
Y <- metadata$Temporada

# 3. Modelo DIABLO Esparso (sPLS-DA)
design <- matrix(0.1, ncol = length(X), nrow = length(X), dimnames = list(names(X), names(X)))
diag(design) <- 0

list.keepX <- list(mRNA = c(30, 30), Proteinas = c(30, 30))
modelo_esparso <- block.splsda(X = X, Y = Y, ncomp = 2, keepX = list.keepX, design = design)

# 4. Visualizaciones Multiómicas
# Circos Plot
png("results/plots/CircosPlot_Integrado.png", width = 1000, height = 1000, res = 120)
circosPlot(modelo_esparso, cutoff = 0.7, line = TRUE, 
           color.blocks = c('darkorchid', 'lightgreen'),
           color.cor = c("chocolate3", "grey20"), size.variables = 0.6)
dev.off()

# Heatmap Multiómico
png("results/plots/Heatmap_Multiomico_CIM.png", width = 1200, height = 1000, res = 100)
cimDiablo(modelo_esparso, margin = c(10, 20), legend.position = "right",
          title = "Heatmap Integrado: Transcriptómica y Proteómica")
dev.off()

# Extracción de Biomarcadores
biomarcadores_arn <- selectVar(modelo_esparso, comp = 1, block = "mRNA")$mRNA$value
biomarcadores_prot <- selectVar(modelo_esparso, comp = 1, block = "Proteinas")$Proteinas$value

write.csv(biomarcadores_arn, "results/tables/Biomarcadores_RNA_Elite.csv")
write.csv(biomarcadores_prot, "results/tables/Biomarcadores_Proteinas_Elite.csv")