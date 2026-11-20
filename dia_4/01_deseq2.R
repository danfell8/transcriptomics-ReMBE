# =============================================================================
# 01_deseq2.R - Expresion diferencial con DESeq2 a partir de kallisto
# Uso en RStudio: corre linea por linea (Ctrl+Enter) o por secciones.
# Necesita: metadata_peces.csv (columnas sample_id, Sample, condition, ...)
#           y txi_<nivel>_deseq2.rds (generado con 00_exportar_counts.R)
# =============================================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)

# -------------------------
# Carpeta de trabajo y parametros
# -------------------------
setwd("~/para_RStudio")

# Nivel del analisis: "transcritos" o "genes"
NIVEL      <- "transcritos"

OUT_DIR    <- paste0("resultados_deseq2_", NIVEL)
REFERENCIA <- "control"     # Nivel de 'condition' contra el que se compara
PADJ       <- 0.05
LFC        <- 1
colores    <- c("steelblue", "firebrick")

dir.create(OUT_DIR, showWarnings = FALSE)

# -------------------------
# 1. Metadatos y cuentas
# -------------------------
meta <- read.csv("metadata_peces.csv")
rownames(meta) <- meta$sample_id
meta$condition <- relevel(factor(meta$condition), ref = REFERENCIA)

txi <- readRDS(paste0("txi_", NIVEL, "_deseq2.rds"))
dim(txi$counts)   # filas = transcritos o genes, columnas = muestras

# MUY IMPORTANTE: metadatos en el mismo orden que las columnas de las cuentas
meta <- meta[colnames(txi$counts), ]
stopifnot(identical(rownames(meta), colnames(txi$counts)))
meta[, c("sample_id", "Sample", "condition")]   # Revisar: 1C-4C control, 5H-8H calor

# -------------------------
# 2. Figura 1: tamaño de las librerias
# -------------------------
barplot(colSums(txi$counts) / 1e6, names.arg = meta$Sample,
        col = colores[meta$condition], las = 1,
        ylab = "Millones de lecturas", main = "Tamaño de las librerias")

# -------------------------
# 3. DESeq2
# -------------------------
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ condition)

# Filtro: al menos 10 cuentas en tantas muestras como el grupo mas pequeno
keep <- rowSums(counts(dds) >= 10) >= min(table(meta$condition))
table(keep)
dds  <- dds[keep, ]

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# -------------------------
# 4. Figura 2: antes y despues de la normalizacion
# -------------------------
# La linea punteada es la mediana de todas las muestras:
# despues de normalizar, las medianas de cada muestra deben quedar alineadas con ella
antes   <- log2(counts(dds) + 1)
despues <- log2(counts(dds, normalized = TRUE) + 1)

par(mfrow = c(1, 2))
boxplot(antes, names = meta$Sample, col = colores[meta$condition],
        las = 2, ylab = "log2(cuentas + 1)", main = "Antes")
abline(h = median(antes), lty = 2, lwd = 2)
boxplot(despues, names = meta$Sample, col = colores[meta$condition],
        las = 2, ylab = "log2(cuentas + 1)", main = "Despues")
abline(h = median(despues), lty = 2, lwd = 2)
par(mfrow = c(1, 1))

# -------------------------
# 5. Figura 3: PCA
# -------------------------
plotPCA(vsd, intgroup = "condition") +
  geom_text(aes(label = meta$Sample), vjust = -1, size = 3) +
  scale_color_manual(values = colores) +
  theme_bw()

# -------------------------
# 6. Resultados
# -------------------------
resultsNames(dds)
coef <- resultsNames(dds)[2]          # "condition_calor_vs_control"

res <- results(dds, name = coef, alpha = PADJ)
summary(res)

res_df <- as.data.frame(res)
res_df$transcrito <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, file.path(OUT_DIR, "DE_todos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig <- subset(res_df, padj < PADJ & abs(log2FoldChange) >= LFC)
write.table(sig, file.path(OUT_DIR, "DE_significativos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
nrow(sig)

# -------------------------
# 7. Figura 4: heatmap de los 50 transcritos mas significativos
# -------------------------
# res_df ya esta ordenado por padj; si hay pocos significativos, usa los 50 primeros
top <- head(res_df$transcrito, 50)
pheatmap(assay(vsd)[top, ], scale = "row", show_rownames = FALSE,
         annotation_col = meta["condition"], labels_col = meta$Sample)

# -------------------------
# 8. Figura 5: volcano plot
# -------------------------
res_df$grupo <- "NS"
res_df$grupo[which(res_df$padj < PADJ & res_df$log2FoldChange >=  LFC)] <- "Up"
res_df$grupo[which(res_df$padj < PADJ & res_df$log2FoldChange <= -LFC)] <- "Down"
table(res_df$grupo)

ggplot(subset(res_df, !is.na(padj)), aes(log2FoldChange, -log10(padj), color = grupo)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed") +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed") +
  theme_bw()