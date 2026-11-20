# =============================================================================
# 00_instalar_paquetes.R - Revisa e instala los paquetes del analisis en R
# Uso: Rscript 00_instalar_paquetes.R   (o corre todo en RStudio)
# Solo instala lo que falte; si ya esta todo, no hace nada.
# =============================================================================

paquetes_cran <- c("ggplot2", "pheatmap")
paquetes_bioc <- c("tximport","edgeR", "DESeq2", "apeglm")

# BiocManager es el instalador de los paquetes de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

faltan_cran <- paquetes_cran[!sapply(paquetes_cran, requireNamespace, quietly = TRUE)]
faltan_bioc <- paquetes_bioc[!sapply(paquetes_bioc, requireNamespace, quietly = TRUE)]

if (length(faltan_cran) > 0) {
  install.packages(faltan_cran, repos = "https://cloud.r-project.org")
}
if (length(faltan_bioc) > 0) {
  BiocManager::install(faltan_bioc, update = FALSE, ask = FALSE)
}

# Reporte final
todos <- c(paquetes_cran, paquetes_bioc)
estado <- sapply(todos, requireNamespace, quietly = TRUE)
print(data.frame(paquete = todos, instalado = estado, row.names = NULL))

if (all(estado)) {
  message("Todo listo.")
} else {
  message("Faltan: ", paste(todos[!estado], collapse = ", "))
}