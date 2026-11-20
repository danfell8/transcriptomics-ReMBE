# =============================================================================
# 02_exploracion_lote.R - ¿De donde viene la variacion entre muestras?
# Correr despues de 01_edgeR.R (usa meta, txi, y, logcpm, colores, OUT_DIR)
#   1. PCA: ¿las muestras se separan por condicion?
#   2. Metadatos: ¿hay alguna otra variable que pueda explicar la dispersion?
#   3. Tamano de libreria: la unica variable tecnica que tenemos a mano
#   4. Un lote hipotetico y tres formas de tratarlo:
#      removeBatchEffect (visualizar), ComBat-seq (corregir cuentas)
#      y el lote dentro del modelo (lo que se haria en un analisis real)
# =============================================================================

library(edgeR)
library(sva)

# Funcion para graficar las dos primeras componentes
graficar_pca <- function(pca, titulo) {
  var_pc <- round(100 * summary(pca)$importance[2, 1:2], 1)
  plot(pca$x[, 1:2], pch = 19, cex = 2, col = colores[meta$condition], main = titulo,
       xlab = paste0("PC1 (", var_pc[1], "%)"),
       ylab = paste0("PC2 (", var_pc[2], "%)"))
  text(pca$x[, 1:2], labels = meta$Sample, pos = 3, cex = 0.7)
}

# -------------------------
# 1. PCA con log2 CPM (calculado en 01_edgeR.R)
# -------------------------
pca <- prcomp(t(logcpm))
graficar_pca(pca, "PCA - log2 CPM")
legend("topright", legend = levels(meta$condition),
       col = colores[1:2], pch = 19)

# Preguntas: ¿PC1 separa control de calor? ¿Hay alguna muestra lejos de su grupo?

# -------------------------
# 2. ¿Que otras variables tenemos en los metadatos?
# -------------------------
# Numero de valores distintos por columna: si es 1, esa variable no puede
# explicar ninguna diferencia entre muestras
sapply(meta, function(x) length(unique(x)))

# -------------------------
# 3. Tamano de libreria (millones de lecturas asignadas por kallisto)
# -------------------------
meta$lib_millones <- colSums(txi$counts) / 1e6
barplot(meta$lib_millones, names.arg = meta$Sample, las = 2,
        col = colores[meta$condition], ylab = "Millones de lecturas")

# ¿Las componentes principales siguen al tamano de libreria?
cor(pca$x[, 1], meta$lib_millones)
cor(pca$x[, 2], meta$lib_millones)

# -------------------------
# 4. Lote hipotetico: las 3 librerias mas profundas
# -------------------------
meta$lote <- factor(ifelse(meta$sample_id %in% c("SRR24971448", "SRR24971449", "SRR24971450"),
                           "B", "A"))
table(meta$lote, meta$condition)   # ¡Ojo! El lote B (6H, 7H, 8H; las mas profundas) es solo calor

# -------------------------
# 5. removeBatchEffect (limma): SOLO para visualizar
# -------------------------
# Quita el efecto del lote de los log2 CPM, protegiendo la condicion.
# No se usa como entrada del analisis estadistico: solo para ver como
# se verian los datos sin el lote.
logcpm_sin_lote <- limma::removeBatchEffect(logcpm, batch = meta$lote,
                                            design = model.matrix(~ condition, data = meta))
pca_l <- prcomp(t(logcpm_sin_lote))

# -------------------------
# 6. ComBat-seq: corrige las CUENTAS
# -------------------------
counts_corr <- ComBat_seq(round(txi$counts), batch = meta$lote, group = meta$condition)

y_c <- DGEList(counts = counts_corr, group = meta$condition)
y_c <- y_c[filterByExpr(y_c), , keep.lib.sizes = FALSE]
y_c <- calcNormFactors(y_c)
pca_c <- prcomp(t(cpm(y_c, log = TRUE)))

# -------------------------
# 7. Comparar los tres PCA lado a lado
# -------------------------
par(mfrow = c(1, 3))
graficar_pca(pca,   "Sin correccion")
graficar_pca(pca_l, "removeBatchEffect")
graficar_pca(pca_c, "ComBat-seq")
par(mfrow = c(1, 1))

# Preguntas: ¿cambia mucho la posicion de las muestras? ¿PC1 explica mas o menos?
# Si el PCA casi no cambia, el lote no es una fuente importante de variacion.

# -------------------------
# 8. Lo que SI se haria en un analisis real: el lote en el modelo
# -------------------------
# No se modifican los datos: se le dice al modelo que hay otra fuente de variacion.
design_lote <- model.matrix(~ lote + condition, data = meta)
colnames(design_lote)   # "(Intercept)" "loteB" "conditioncalor"

y_l   <- estimateDisp(y, design_lote)
qlf_l <- glmQLFTest(glmQLFit(y_l, design_lote), coef = "conditioncalor")

# ComBat-seq en el modelo de siempre (~ condition)
design <- model.matrix(~ condition, data = meta)
y_c    <- estimateDisp(y_c, design)
qlf_c  <- glmQLFTest(glmQLFit(y_c, design), coef = 2)

# ¿Cambia el numero de transcritos diferenciales?
summary(decideTests(qlf))     # Sin correccion (de 01_edgeR.R)
summary(decideTests(qlf_l))   # Con el lote en el modelo
summary(decideTests(qlf_c))   # Con ComBat-seq
