# =============================================================================
# 00_exportar_counts.R - Junta las salidas de kallisto en objetos de R
# Se corre en el servidor del curso (donde estan las carpetas de kallisto):
#   sbatch R_analisis.sh
# Luego se copian a RStudio solo los archivos de la carpeta para_RStudio/
# =============================================================================

library(tximport)

# dropInfReps = TRUE: no lee los bootstraps del abundance.h5 (evita necesitar rhdf5)

KALLISTO_DIR <- "kallisto_conteos"
METADATOS    <- "../metadata_peces.csv"
OUT_DIR      <- "para_RStudio"

dir.create(OUT_DIR, showWarnings = FALSE)

meta  <- read.csv(METADATOS)
files <- file.path(KALLISTO_DIR, meta$sample_id, "abundance.tsv")
names(files) <- meta$sample_id
# Si falta algun archivo, dice cual (revisa KALLISTO_DIR y los sample_id)
faltan <- files[!file.exists(files)]
if (length(faltan) > 0) stop("No encuentro estos archivos:\n", paste(faltan, collapse = "\n"))

# Para DESeq2: cuentas estimadas + longitudes de los transcritos
txi <- tximport(files, type = "kallisto", txOut = TRUE, dropInfReps = TRUE)
saveRDS(txi, file.path(OUT_DIR, "txi_deseq2.rds"))

# Para edgeR: cuentas corregidas por longitud
txi_edger <- tximport(files, type = "kallisto", txOut = TRUE, dropInfReps = TRUE,
                      countsFromAbundance = "lengthScaledTPM")
saveRDS(txi_edger, file.path(OUT_DIR, "txi_edgeR.rds"))

# Tablas en texto, por si se quieren revisar en Excel o en otro programa
write.table(data.frame(transcrito = rownames(txi$counts), txi$counts, check.names = FALSE),
            file.path(OUT_DIR, "counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(transcrito = rownames(txi$abundance), txi$abundance, check.names = FALSE),
            file.path(OUT_DIR, "TPM.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

file.copy(METADATOS, OUT_DIR, overwrite = TRUE)

message("Listo. Copia la carpeta ", OUT_DIR, "/ al servidor de RStudio.")