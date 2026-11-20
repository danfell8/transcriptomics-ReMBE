# =============================================================================
# 01_edgeR.R - Expresion diferencial con edgeR a partir de kallisto
# Uso en RStudio: corre linea por linea (Ctrl+Enter) o por secciones.
# Necesita: metadata_peces.csv (columnas sample_id, Sample, condition, ...)
#           y txi_<nivel>_edgeR.rds (generado con 00_exportar_counts.R)
# =============================================================================

library(edgeR)
library(ggplot2)
library(pheatmap)

# -------------------------
# Carpeta de trabajo y parametros
# -------------------------
setwd("~/para_RStudio")

# Nivel del analisis: "transcritos" o "genes"
NIVEL      <- "transcritos"

OUT_DIR    <- paste0("resultados_edgeR_", NIVEL)
REFERENCIA <- "control"     # Nivel de 'condition' contra el que se compara
UMBRAL_FDR <- 0.05
LFC        <- 1
colores    <- c("steelblue", "firebrick")

dir.create(OUT_DIR, showWarnings = FALSE)

# -------------------------
# 1. Metadatos y cuentas
# -------------------------
meta <- read.csv("metadata_peces.csv")
rownames(meta) <- meta$sample_id
meta$condition <- relevel(factor(meta$condition), ref = REFERENCIA)

txi <- readRDS(paste0("txi_", NIVEL, "_edgeR.rds"))
dim(txi$counts)   # filas = transcritos o genes, columnas = muestras

# MUY IMPORTANTE: los metadatos deben estar en el mismo orden que las columnas
# de la matriz de cuentas. Si no, cada muestra queda con la etiqueta de otra.
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
# 3. edgeR
# -------------------------
y <- DGEList(counts = txi$counts, group = meta$condition)

# Filtro de transcritos con muy pocas cuentas
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalizacion TMM
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE)

# -------------------------
# 4. Figura 2: antes y despues de la normalizacion
# -------------------------
# La linea punteada es la mediana de todas las muestras:
# despues de normalizar, las medianas de cada muestra deben quedar alineadas con ella
antes <- log2(y$counts + 1)

par(mfrow = c(1, 2))
boxplot(antes, names = meta$Sample, col = colores[meta$condition],
        las = 2, ylab = "log2(cuentas + 1)", main = "Antes")
abline(h = median(antes), lty = 2, lwd = 2)
boxplot(logcpm, names = meta$Sample, col = colores[meta$condition],
        las = 2, ylab = "log2 CPM (TMM)", main = "Despues")
abline(h = median(logcpm), lty = 2, lwd = 2)
par(mfrow = c(1, 1))

# -------------------------
# 5. Figura 3: PCA
# -------------------------
pca    <- prcomp(t(logcpm))
var_pc <- round(100 * summary(pca)$importance[2, 1:2], 1)
plot(pca$x[, 1:2], pch = 19, cex = 2, col = colores[meta$condition],
     xlab = paste0("PC1 (", var_pc[1], "%)"),
     ylab = paste0("PC2 (", var_pc[2], "%)"), main = "PCA")
text(pca$x[, 1:2], labels = meta$Sample, pos = 3)
legend("topright", legend = levels(meta$condition), col = colores, pch = 19)

# -------------------------
# 6. Modelo y resultados
# -------------------------
design <- model.matrix(~ condition, data = meta)
y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

colnames(design)
qlf <- glmQLFTest(fit, coef = 2)       # Columna 2: condicion vs referencia
summary(decideTests(qlf))

res_df <- as.data.frame(topTags(qlf, n = Inf))
res_df$transcrito <- rownames(res_df)
write.table(res_df, file.path(OUT_DIR, "DE_todos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig <- subset(res_df, FDR < UMBRAL_FDR & abs(logFC) >= LFC)
write.table(sig, file.path(OUT_DIR, "DE_significativos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
nrow(sig)

# -------------------------
# 7. Figura 4: heatmap de los 50 transcritos mas significativos
# -------------------------
# res_df ya esta ordenado por FDR; si hay pocos significativos, usa los 50 primeros
top <- head(res_df$transcrito, 50)
pheatmap(logcpm[top, ], scale = "row", show_rownames = FALSE,
         annotation_col = meta["condition"], labels_col = meta$Sample)

# -------------------------
# 8. Figura 5: volcano plot
# -------------------------
res_df$grupo <- "NS"
res_df$grupo[res_df$FDR < UMBRAL_FDR & res_df$logFC >=  LFC] <- "Up"
res_df$grupo[res_df$FDR < UMBRAL_FDR & res_df$logFC <= -LFC] <- "Down"
table(res_df$grupo)

ggplot(res_df, aes(logFC, -log10(FDR), color = grupo)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed") +
  geom_hline(yintercept = -log10(UMBRAL_FDR), linetype = "dashed") +
  theme_bw()