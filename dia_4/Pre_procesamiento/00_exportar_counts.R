# =============================================================================
# 00_exportar_counts.R - Junta las salidas de kallisto en objetos de R
# Se corre en el servidor del curso:  sbatch R_analisis.sh
# Exporta las cuentas a DOS niveles:
#   - transcritos: una fila por transcrito (como los cuantifico kallisto)
#   - genes:       se suman los transcritos (isoformas) de cada gen
# Luego se copian a RStudio solo los archivos de la carpeta para_RStudio/
# =============================================================================

library(tximport)

KALLISTO_DIR <- "kallisto_conteos"
METADATOS    <- "../metadata_peces.csv"
FASTA        <- "/scratch/taller/datos_pez/transcriptoma_completo/GKNI01.1.fsa_nt.gz"   # el mismo que se uso para kallisto
OUT_DIR      <- "para_RStudio"

dir.create(OUT_DIR, showWarnings = FALSE)

# -------------------------
# 1. Archivos de kallisto
# -------------------------
meta  <- read.csv(METADATOS)
files <- file.path(KALLISTO_DIR, meta$sample_id, "abundance.tsv")
names(files) <- meta$sample_id

# Si falta algun archivo, dice cual (revisa KALLISTO_DIR y los sample_id)
faltan <- files[!file.exists(files)]
if (length(faltan) > 0) stop("No encuentro estos archivos:\n", paste(faltan, collapse = "\n"))

# -------------------------
# 2. Nivel TRANSCRITO (txOut = TRUE: no se agrupa nada)
# -------------------------
# dropInfReps = TRUE: no lee los bootstraps del abundance.h5 (evita necesitar rhdf5)

# Para DESeq2: cuentas estimadas + longitudes de los transcritos
txi_tx <- tximport(files, type = "kallisto", txOut = TRUE, dropInfReps = TRUE)
saveRDS(txi_tx, file.path(OUT_DIR, "txi_transcritos_deseq2.rds"))

# Para edgeR: cuentas corregidas por longitud
txi_tx_edger <- tximport(files, type = "kallisto", txOut = TRUE, dropInfReps = TRUE,
                         countsFromAbundance = "lengthScaledTPM")
saveRDS(txi_tx_edger, file.path(OUT_DIR, "txi_transcritos_edgeR.rds"))

# -------------------------
# 3. Tabla transcrito -> gen (desde los encabezados del FASTA)
# -------------------------
# Encabezado: >GKNI01000001.1 TSA: Genypterus maculatus CONTIG_10000_c0_g2_i1.p1, ...
#   transcrito = GKNI01000001.1   (el ID que uso kallisto)
#   gen        = CONTIG_10000_c0_g2   (sin la isoforma _i1 ni el ORF .p1)
encabezados <- sub("^>", "", grep("^>", readLines(gzfile(FASTA)), value = TRUE))
patron      <- "CONTIG_[0-9]+_c[0-9]+_g[0-9]+"

tx2gene <- data.frame(
  transcrito = sub(" .*", "", encabezados),
  gen        = ifelse(grepl(patron, encabezados),
                      sub(paste0(".*(", patron, ").*"), "\\1", encabezados),
                      sub(" .*", "", encabezados))   # sin nombre de gen: queda solo
)

stopifnot(all(rownames(txi_tx$counts) %in% tx2gene$transcrito))
message("Transcritos: ", nrow(tx2gene), "  Genes: ", length(unique(tx2gene$gen)))
write.table(tx2gene, file.path(OUT_DIR, "tx2gene.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 4. Nivel GEN (se suman los transcritos de cada gen)
# -------------------------
txi_gen <- summarizeToGene(txi_tx, tx2gene)
saveRDS(txi_gen, file.path(OUT_DIR, "txi_genes_deseq2.rds"))

txi_gen_edger <- summarizeToGene(txi_tx, tx2gene, countsFromAbundance = "lengthScaledTPM")
saveRDS(txi_gen_edger, file.path(OUT_DIR, "txi_genes_edgeR.rds"))

# -------------------------
# 5. Tablas en texto y metadatos
# -------------------------
write.table(data.frame(transcrito = rownames(txi_tx$counts), txi_tx$counts, check.names = FALSE),
            file.path(OUT_DIR, "counts_transcritos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gen = rownames(txi_gen$counts), txi_gen$counts, check.names = FALSE),
            file.path(OUT_DIR, "counts_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

file.copy(METADATOS, OUT_DIR, overwrite = TRUE)

message("Listo. Copia la carpeta ", OUT_DIR, "/ al servidor de RStudio.")