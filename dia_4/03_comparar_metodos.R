# =============================================================================
# 03_comparar_metodos.R - ¿Cambia el resultado segun el metodo estadistico?
# Script independiente: se puede correr en una sesion limpia.
# Mismos datos, cinco pruebas: edgeR (QL, QL robusta, LRT, exacta) y DESeq2.
# Criterio igual para todas: FDR (padj) < 0.05, sin corte de log2FC.
# =============================================================================

library(edgeR)
library(DESeq2)

setwd("~/para_RStudio")
FDR   <- 0.05
NIVEL <- "transcritos"   # o "genes"

# -------------------------
# 1. Metadatos y cuentas
# -------------------------
meta <- read.csv("metadata_peces.csv")
rownames(meta) <- meta$sample_id
meta$condition <- relevel(factor(meta$condition), ref = "control")

txi_edger  <- readRDS(paste0("txi_", NIVEL, "_edgeR.rds"))
txi_deseq2 <- readRDS(paste0("txi_", NIVEL, "_deseq2.rds"))

meta <- meta[colnames(txi_edger$counts), ]
stopifnot(identical(rownames(meta), colnames(txi_edger$counts)),
          identical(rownames(meta), colnames(txi_deseq2$counts)))

design <- model.matrix(~ condition, data = meta)

# -------------------------
# 2. edgeR: cuatro pruebas sobre el mismo objeto
# -------------------------
y <- DGEList(counts = txi_edger$counts, group = meta$condition)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

sig_edgeR <- function(prueba) rownames(topTags(prueba, n = Inf, p.value = FDR))

sig_ql     <- sig_edgeR(glmQLFTest(glmQLFit(y, design), coef = 2))
sig_ql_rob <- sig_edgeR(glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = 2))
sig_lrt    <- sig_edgeR(glmLRT(glmFit(y, design), coef = 2))
sig_exacta <- sig_edgeR(exactTest(y))

# -------------------------
# 3. DESeq2
# -------------------------
dds  <- DESeqDataSetFromTximport(txi_deseq2, colData = meta, design = ~ condition)
dds  <- dds[rowSums(counts(dds) >= 10) >= min(table(meta$condition)), ]
dds  <- DESeq(dds)
res  <- results(dds, alpha = FDR)
sig_deseq2 <- rownames(subset(as.data.frame(res), padj < FDR))

# -------------------------
# 4. ¿Cuantos transcritos encuentra cada metodo?
# -------------------------
n_sig <- c(QL         = length(sig_ql),
           QL_robusta = length(sig_ql_rob),
           Exacta     = length(sig_exacta),
           LRT        = length(sig_lrt),
           DESeq2     = length(sig_deseq2))
n_sig

barplot(n_sig, col = c(rep("steelblue", 4), "firebrick"), las = 1,
        ylab = "Transcritos diferenciales (FDR < 0.05)",
        main = "Mismos datos, distinto metodo")

# -------------------------
# 5. ¿Encuentran los MISMOS transcritos?
# -------------------------
# Diagrama de Venn (QL no se incluye porque no encontro ningun transcrito)
todos <- unique(c(sig_exacta, sig_lrt, sig_deseq2))
venn  <- cbind(Exacta = todos %in% sig_exacta,
               LRT    = todos %in% sig_lrt,
               DESeq2 = todos %in% sig_deseq2)
limma::vennDiagram(limma::vennCounts(venn),
                   circle.col = c("steelblue", "darkgreen", "firebrick"))

nucleo <- Reduce(intersect, list(sig_exacta, sig_lrt, sig_deseq2))
length(nucleo)   # Transcritos que encuentran los tres metodos

# Tabla del nucleo con los resultados de DESeq2
as.data.frame(res[nucleo, c("baseMean", "log2FoldChange", "padj")])

# Preguntas para discutir:
#  - ¿Por que QL no encuentra nada si los otros metodos si?
#  - ¿Que transcritos reportarias? ¿Por que?
#  - ¿Cuando se debe elegir el metodo: antes o despues de ver los resultados?