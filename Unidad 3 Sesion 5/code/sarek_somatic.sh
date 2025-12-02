#!/bin/bash
# Ejecuta nf-core/sarek en modo SOMÁTICO (tumor-only) para la muestra S9.
# Basado en el script guía entregado por la profesora.
# El script crea internamente un samplesheet CSV como requiere nf-core/sarek.

###############################
# Definir entradas para S9    #
###############################

# Rutas a los FASTQ crudos de la muestra S9
R1="/home/bioinfo1/181004_curso_calidad_datos_NGS/fastq_raw/S9_R1.fastq.gz"
R2="/home/bioinfo1/181004_curso_calidad_datos_NGS/fastq_raw/S9_R2.fastq.gz"

# Directorio de salida (desde la carpeta code/)
OUT="../results/somatic_S9"

# Crear el directorio de salida
mkdir -p "$OUT"

############################################
# Detección automática del nombre (S9)     #
############################################

base=$(basename "$R1")

# Elimina sufijos comunes de R1 (igual que en el script original)
sample=${base%%_R1.fastq.gz}
sample=${sample%%_R1.fq.gz}
sample=${sample%%_1.fastq.gz}
sample=${sample%%_1.fq.gz}
sample=${sample%%.fastq.gz}
sample=${sample%%.fq.gz}

SAMPLE=$sample
echo "Detectado nombre de muestra automáticamente: ${SAMPLE}"

############################################
# Rutas absolutas                          #
############################################

if command -v readlink >/dev/null 2>&1; then
    R1_ABS=$(readlink -f "$R1")
    R2_ABS=$(readlink -f "$R2")
else
    R1_ABS="$R1"
    R2_ABS="$R2"
fi

############################################
# Crear el samplesheet CSV para Sarek      #
############################################

SHEET="${OUT}/samplesheet_somatic_${SAMPLE}.csv"

echo "Creando samplesheet: $SHEET"
cat > "$SHEET" <<EOF
patient,sex,status,sample,lane,fastq_1,fastq_2
${SAMPLE},NA,1,${SAMPLE},L1,${R1_ABS},${R2_ABS}
EOF

############################################
# Ejecutar nf-core/sarek en modo somático  #
############################################

echo "Lanzando nf-core/sarek en modo somático (tumor-only) para ${SAMPLE}..."
nextflow run nf-core/sarek \
    --input "$SHEET" \
    --genome GATK.GRCh38 \
    --outdir "$OUT" \
    --tools mutect2 \
    -profile singularity \
    -c /home/bioinfo1/korostica/test_tutorial/code/local_sarek_8cpus.config \
    -resume

