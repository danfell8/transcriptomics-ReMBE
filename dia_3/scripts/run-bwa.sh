#!/bin/bash
#SBATCH -J bwa-oyamel                 # Nombre del trabajo
#SBATCH -p short                      # Partición (short para pruebas/taller)
#SBATCH -c 4                          # Número de hilos (CPUs por tarea)
#SBATCH --mem=16G                     # Memoria RAM solicitada (16G para sort y BWA)
#SBATCH -o ../error-data-outs/bwa-oyamel_%j.out      # Archivo de salida estándar
#SBATCH -e ../error-data-outs/bwa-oyamel_%j.err      # Archivo de errores

# ==============================================================================
# PASO 1: Iniciar Conda y activar el ambiente compartido
# ==============================================================================

echo ">>> Cargando módulos y activando entorno rnaseq_env..."
module load spack/1.0.2
module load miniconda3/24.7.1-none-none-linux-almalinux9-x86_64-g6mq4vw

eval "$(conda shell.bash hook)"
conda activate /home/veronica/.conda/envs/rnaseq_env

# ==============================================================================
# 2. Definir rutar y archivos de salida
# ==============================================================================
# Datos compartidos (solo lectura)

REF_FASTA="/scratch/taller/datos_oyamel/REFERENCE/GCAT_AB-RNA-1.0.16.fa"
TRIMMED_DIR="/scratch/taller/quality-assessment/results/trimmed"

# Directorios de salida en la carpeta local 'outputs/' de cada usuari@

OUT_BASE="../outputs/BWA"
INDEX_DIR="${OUT_BASE}/Index_bwa"
INDEX_PREFIX="${INDEX_DIR}/index_Areligiosa"
SAM_DIR="${OUT_BASE}/SAM"
BAM_DIR="${OUT_BASE}/BAM"
COUNT_DIR="${OUT_BASE}/TXT/COUNT_GENES"
ORDER_DIR="${OUT_BASE}/TXT/GENES_ORDER"

# Crear carpetas

mkdir -p "${INDEX_DIR}" "${SAM_DIR}" "${BAM_DIR}" "${COUNT_DIR}" "${ORDER_DIR}"

# ==============================================================================
# 3. Construir índice de BWA
# ==============================================================================

echo ">>> [Paso 1/4] Construyendo índice de BWA para ${REF_FASTA}..."
bwa index -p "${INDEX_PREFIX}" "${REF_FASTA}"

# ==============================================================================
# 4. Alineamiento, conversión a BAM y conteo por muestra
# ==============================================================================
# Lista de muestras del taller

MUESTRAS=("DPVR1_S179" "DPVR3_S181" "DPVR4_S182" "DPVR6_S184" "DPVR8_S186" "DPVR9_S187")

for a in "${MUESTRAS[@]}"; do
    echo "=========================================================="
    echo "Procesando muestra: ${a}"
    echo "=========================================================="

    # Rutas a los archivos pareados limpios 
    R1="${TRIMMED_DIR}/Trimmer_${a}_L007_R1_001_paired.fq.gz"
    R2="${TRIMMED_DIR}/Trimmer_${a}_L007_R2_001_paired.fq.gz"

    # Archivos de salida
    SAM_FILE="${SAM_DIR}/${a}_sw10L50.sam"
    BAM_FILE="${BAM_DIR}/${a}_sw10L50.bam"
    COUNT_FILE="${COUNT_DIR}/${a}_sw10L50.countgenes.txt"
    ORDER_FILE="${ORDER_DIR}/${a}_sw10L50.genesorder.txt"

    # --- PASO 1: Mapeo con BWA MEM ---
    echo ">>> [Paso 1/3] Alineando lecturas pareadas con BWA MEM..."
    bwa mem -t "${SLURM_CPUS_PER_TASK}" \
        "${INDEX_PREFIX}" \
        "${R1}" \
        "${R2}" \
        > "${SAM_FILE}"

    # --- PASO 2: Conversión de SAM a BAM ---
    echo ">>> [Paso 2/3] Convirtiendo SAM a BAM comprimido..."
    samtools view -@ "${SLURM_CPUS_PER_TASK}" -bS "${SAM_FILE}" > "${BAM_FILE}"

    # (Opcional, pero altamente recomendable) Eliminar el archivo .sam pesado una vez creado el .bam para ahorrar espacio:
    # rm -f "${SAM_FILE}"

    # --- PASO 3: Conteo de lecturas mapeadas por transcrito ---
    echo ">>> [Paso 3/3] Extrayendo y contando transcritos mapeados (-F 260)..."
    samtools view -@ "${SLURM_CPUS_PER_TASK}" -F 260 "${BAM_FILE}" \
        | cut -f 3 \
        | sort --parallel="${SLURM_CPUS_PER_TASK}" -S 8G \
        | uniq -c \
        | tee "${COUNT_FILE}" \
        | awk '{print $1 "\t" $2}' > "${ORDER_FILE}"

    echo ">>> Muestra ${a} finalizada con éxito."
done

echo "=========================================================="
echo ">>> ¡Pipeline de BWA y conteo completado para todas las muestras!"
