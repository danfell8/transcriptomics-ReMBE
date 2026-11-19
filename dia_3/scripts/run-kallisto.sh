#!/bin/bash
#SBATCH -J kallisto-oyamel            # Nombre del trabajo
#SBATCH -p short                      # Partición (short para pruebas/taller)
#SBATCH -c 4                          # Número de hilos (CPUs por tarea)
#SBATCH --mem=8G                     # Memoria RAM solicitada (16G para sort y BWA)
#SBATCH -o ../error-data-outs/kallisto-oyamel_%j.out      # Archivo de salida estándar
#SBATCH -e ../error-data-outs/kallisto-oyamel_%j.err      # Archivo de errores

# ==============================================================================
# PASO 1: Iniciar Conda y activar el ambiente compartido
# ==============================================================================

echo ">>> Cargando módulos y activando entorno rnaseq_env..."
module load spack/1.0.2
module load miniconda3/24.7.1-none-none-linux-almalinux9-x86_64-g6mq4vw

eval "$(conda shell.bash hook)"
conda activate /home/veronica/.conda/envs/rnaseq_env

# ==============================================================================
# PASO 2: Definir rutas de entrada y salida
# ==============================================================================
# El script se guarda y ejecuta desde la carpeta scripts/ que cada usuari@ clono con GitHub

# Datos compartidos (solo lectura)

REF_FASTA="/scratch/taller/datos_oyamel/REFERENCE/GCAT_AB-RNA-1.0.16.fa"
TRIM_DIR="/scratch/taller/quality-assessment/results/trimmed"

# Directorios de salida en la carpeta local 'outputs/' de cada usuari@

OUT_DIR="../outputs/KALLISTO"
INDEX_FILE="${OUT_DIR}/GCAT_AB-RNA-1.0.16.idx"

mkdir -p "${OUT_DIR}"

# ==============================================================================
# PASO 3: Construir el índice del transcriptoma de referencia
# ==============================================================================

echo "Construyendo índice de Kallisto..."
kallisto index -i "${INDEX_FILE}" "${REF_FASTA}"

# ==============================================================================
# PASO 4: Cuantificar cada muestra paired-end
# ==============================================================================
for R1 in "${TRIM_DIR}"/*_R1_001_paired.fq.gz; do
    # Extraer el prefijo de la muestra (ej. Trimmer_DPVR1_S179_L007)
    SAMPLE=$(basename "${R1}" _R1_001_paired.fq.gz)
    R2="${TRIM_DIR}/${SAMPLE}_R2_001_paired.fq.gz"

    echo "Cuantificando muestra: ${SAMPLE}... finalizada con éxito."
    kallisto quant \
        -i "${INDEX_FILE}" \
        -o "${OUT_DIR}/${SAMPLE}" \
        -b 100 \
        -t 4 \
	"${R1}" "${R2}"
done

echo "=========================================================="
echo ">>>¡Pipeline de Kallisto y conteo completado para todas las muestras!"
