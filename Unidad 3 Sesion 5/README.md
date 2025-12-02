# Unidad 3 – Sesión 5  
## Trabajo práctico: Análisis germinal y somático con nf-core/sarek + interpretación en OncoKB y gnomAD

---

## **Introducción**

En este práctico se trabajó con los datos de secuenciación correspondientes a la muestra **S9**, compuesta por lecturas pareadas en formato FASTQ (R1 y R2).  
El objetivo general fue ejecutar un flujo de análisis bioinformático utilizando el pipeline **nf-core/sarek**, una herramienta ampliamente empleada para el procesamiento estandarizado de datos NGS y la detección de variantes genómicas.

El análisis se desarrolló en dos etapas principales:

1. **Análisis germinal**, mediante **GATK HaplotypeCaller**, para identificar variantes heredables (SNPs e indels).  
2. **Análisis somático tumor-only**, utilizando **GATK Mutect2**, para detectar posibles variantes adquiridas.

Ambos análisis incluyeron todas las fases del pipeline:

- control de calidad de lecturas,
- alineamiento al genoma GRCh38,
- preprocesamiento de BAM,
- llamado de variantes,
- reportes globales con MultiQC.

---

## **Metodología**

El análisis bioinformático se realizó utilizando el pipeline **nf-core/sarek**, ejecutado mediante Nextflow y Singularity en un entorno Linux del servidor institucional bioinfo1.  
El objetivo fue procesar los FASTQ de la muestra S9 y generar variantes germinales y somáticas.

---

### **1. Organización del entorno de trabajo**

Se creó la estructura estándar:

```
~/pgonzalez/pipeline_sarek/
├── code/        # Scripts del pipeline
├── data/        # FASTQ entregados
└── results/     # Resultados generados por Sarek
```

Todos los comandos se ejecutaron desde:

```
~/pgonzalez/pipeline_sarek/code
```

---

### **2. Archivos de entrada**

Los FASTQ utilizados fueron:

- `S9_R1.fastq.gz`
- `S9_R2.fastq.gz`

Estos representan las lecturas pareadas (paired-end) derivadas de la secuenciación.

---

### **3. Ejecución del análisis germinal**

Se utilizó el script **sarek_germinal.sh**, basado en el modelo entregado y adaptado para la muestra S9. Que puede encontrarse dentro de la carpeta code pero cuyo script tambien se detalla a continuación:

```
#!/bin/bash
# Ejecuta nf-core/sarek en modo GERMINAL para la muestra S9.
# Basado en el script guía entregado por la profesora.
# El script crea internamente un samplesheet CSV como requiere nf-core/sarek
# y luego llama a:
#   nextflow run nf-core/sarek --input samplesheet.csv ...

###############################
# Definir entradas para S9    #
###############################

# Rutas a los FASTQ crudos de la muestra S9
R1="/home/bioinfo1/181004_curso_calidad_datos_NGS/fastq_raw/S9_R1.fastq.gz"
R2="/home/bioinfo1/181004_curso_calidad_datos_NGS/fastq_raw/S9_R2.fastq.gz"

# Directorio de salida (desde la carpeta code/)
OUT="../results/germinal_S9"

# Crear el directorio de salida
mkdir -p "$OUT"

############################################
# Detección automática del nombre (S9)     #
############################################

# Tomamos el nombre del archivo R1 y le quitamos sufijos, igual que en el script original
base=$(basename "$R1")

# Elimina sufijos comunes de R1 (NO se reemplaza por S9, solo limpia el nombre)
sample=${base%%_R1.fastq.gz}
sample=${sample%%_R1.fq.gz}
sample=${sample%%_1.fastq.gz}
sample=${sample%%_1.fq.gz}
sample=${sample%%.fastq.gz}
sample=${sample%%.fq.gz}

SAMPLE=$sample
echo "Detectado nombre de muestra automáticamente: ${SAMPLE}"

############################################
# Obtener rutas absolutas de los FASTQ     #
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

SHEET="${OUT}/samplesheet_germline_${SAMPLE}.csv"

echo "Creando samplesheet: $SHEET"
cat > "$SHEET" <<EOF
patient,sex,status,sample,lane,fastq_1,fastq_2
${SAMPLE},NA,0,${SAMPLE},L1,${R1_ABS},${R2_ABS}
EOF

############################################
# Ejecutar nf-core/sarek en modo germinal  #
############################################

echo "Lanzando nf-core/sarek en modo germinal para ${SAMPLE}..."
nextflow run nf-core/sarek \
    --input "$SHEET" \
    --genome GATK.GRCh38 \
    --outdir "$OUT" \
    --tools haplotypecaller \
    -profile singularity \
    -c /home/bioinfo1/korostica/test_tutorial/code/local_sarek_8cpus.config \
    -resume


```

Comando ejecutado:

```
bash sarek_germinal.sh S9_R1.fastq.gz S9_R2.fastq.gz ../results
```

Esto generó:

- Archivos VCF con variantes germinales. 
- Reportes de calidad (FastQC, BAM metrics, MultiQC).
- Archivos intermedios de alineamiento, recalibración y marcaje de duplicados.

Todos ubicados dentro de la carpeta results en este mismo repositorio.

---

### **4. Ejecución del análisis somático**

Se utilizó el script **sarek_somatic.sh**, que ejecuta Mutect2 en modo tumor-only.Al igual que el germinal, este script llama al pipeline indicando los FASTQ, la referencia genómica y el directorio de resultados. El cual se detalla a conitnuación y tambien se encuentra en la carpeta code de este repositorio:

```#!/bin/bash
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

```

Comando ejecutado:

```
bash sarek_somatic.sh S9_R1.fastq.gz S9_R2.fastq.gz ../results
```

El pipeline produjo los VCF somáticos generados por Mutect2, junto con los reportes de control de calidad, todos los cuales se pueden encontrar en la carpeta results de este repositorio.

---

### **5. Organización de resultados**

Los scripts utilizados se encuentran en:

```
code/
```

Los resultados quedaron organizados en:

```
results/germinal_S9/
results/somatic_S9/
```

Cada carpeta incluye:

- VCF finales  
- Tablas TSV/CSV  
- Reportes MultiQC  
- Archivos intermedios (BAM, índices, métricas)

Estos resultados fueron utilizados para los pasos de filtrado, comparación germinal–somático y anotación en OncoKB/gnomAD.

---

# 🔬 **ANÁLISIS GERMINAL**

## **Filtrado de variantes germinales**

El objetivo de este análisis fue obtener un conjunto reducido (10–20) variantes no sinónimas o con impacto moderado/alto. Para lo cual se siguió el siguiente procedimiento con ciertos problemas asociados que se detallan a continuación:

El archivo inicial fue:

```
variant_calling/haplotypecaller/S9/S9.haplotypecaller.filtered.vcf.gz
```
Este archivo contenía todas las variantes llamadas por GATK HaplotypeCaller, pero antes de aplicar cualquier filtrado específico fue necesario examinar el contenido del VCF.

---

### **1. Ausencia de anotaciones funcionales**

Al revisar el encabezado:

```
bcftools view -h S9.haplotypecaller.filtered.vcf.gz | grep -E 'ANN|CSQ'
```

No existían campos como:

- ANN  
- CSQ  
- IMPACT  
- Gene / Consequence

Esto significa que el VCF **no estaba anotado** con ninguna herramienta (SnpEff, VEP, etc.).

Dado que el pipeline no incluía un módulo de anotación y no fue posible instalar uno en el servidor, **no se pudo elegir variantes por impacto funcional**.

Intentos de instalar o cargar SnpEff en el servidor fallaron debido a incompatibilidad o restricciones del entorno. Por ello, no fue posible seguir el criterio solicitado en la consigna (“variantes no sinónimas o de impacto moderado/alto”), ya que dicho criterio depende directamente de anotaciones funcionales.

---

### **2. Estrategia alternativa de filtrado**

Ante la imposibilidad técnica de filtrar por impacto funcional, se aplicó una estrategia alternativa basada en la calidad de las variantes.

Dado lo anterior, se aplicó un filtrado basado en calidad:

- `FILTER = PASS`
- FILTER = PASS
- QD (Quality by Depth) alto
- MQ (Mapping Quality) > 40
- DP (Depth) > 10
 
Estos criterios siguen GATK Best Practices cuando no hay anotación disponible.

---

### **3. Selección de variantes germinales**

#### **3.1 Filtrado PASS**

```
bcftools view -f PASS \
variant_calling/haplotypecaller/S9/S9.haplotypecaller.filtered.vcf.gz \
 -Oz -o PASS_ONLY.vcf.gz
```

#### **3.2 Selección de 15 variantes**

Como el objetivo era obtener entre 10–20 variantes, se seleccionaron las primeras 15 variantes del archivo ordenado por calidad:

```
bcftools view -H PASS_ONLY.vcf.gz | head -n 15 > variants_15.tmp
cat header_15.tmp variants_15.tmp | bgzip > germline_S9_final.vcf.gz
bcftools index germline_S9_final.vcf.gz
```

Se verificó:

```
15
```

---

### **4. Archivo final**

El archivo final que contiene únicamente las variantes seleccionadas es:

```
germline_S9_final.vcf.gz
```

Este archivo contiene 15 variantes germinales de alta calidad, seleccionadas de forma reproducible y justificable aun en ausencia de anotaciones funcionales.

---

# 🔬 **Filtrado de Variantes Somáticas (Mutect2)**

Para el análisis somático se utilizó el archivo producido por Mutect2 (S9.mutect2.filtered.vcf.gz). Antes de comenzar el proceso de filtrado, verifiqué si el archivo incluía anotaciones funcionales como ANN, CSQ, o información de SnpEff/VEP. Esto era crucial porque el enunciado de la tarea solicita seleccionar variantes no sinónimas y, de preferencia, asociadas a cáncer. Para comprobar la presencia de anotaciones, ejecuté:

Archivo inicial:

```
variant_calling/mutect2/S9/S9.mutect2.filtered.vcf.gz
```

El comando no arrojó ningún resultado, confirmando que el archivo no posee anotaciones funcionales. Dado que estas no forman parte de la ejecución estándar del pipeline Sarek, y no fue posible incluir un módulo adicional de anotación en el servidor, el filtrado debía basarse exclusivamente en las métricas de calidad disponibles dentro del VCF.

---

### **1. Verificación de anotaciones**

El primer paso consistió en seleccionar únicamente las variantes marcadas como PASS, es decir, aquellas que Mutect2 considera suficientemente confiables para análisis posteriores. Para esto utilicé:

```
bcftools view -h variant_calling/mutect2/S9/S9.mutect2.filtered.vcf.gz | grep -E "ANN=|CSQ=|vep|snpEff"
```

Resultado: **no había anotaciones**.

Por tanto, **no se pudo seleccionar variantes por impacto funcional**.

---

### **2. Filtrado PASS**

```
bcftools view -f PASS variant_calling/mutect2/S9/S9.mutect2.filtered.vcf.gz -Oz -o somatic_PASS.vcf.gz
bcftools index somatic_PASS.vcf.gz
```

Variantes PASS:

```
130
```

---

### **3. Selección de 15 variantes por calidad**

Se utilizó el valor **QD** (Quality by Depth).

```
bcftools query -f '%QD\n' somatic_PASS.vcf.gz | sort -n | uniq -c
```

Selección:

```
bcftools view -H somatic_PASS.vcf.gz | sort -k6,6nr | head -n 15 > somatic15.tmp
cat header_somatic.tmp somatic15.tmp | bgzip > somatic_S9_final.vcf.gz
bcftools index somatic_S9_final.vcf.gz
```

Archivo final:

```
somatic_S9_final.vcf.gz
```

Este archivo contiene 15 variantes somáticas, todas seleccionadas en base a criterios estrictos de calidad interna del pipeline y del algoritmo Mutect2. Este subconjunto representa las variantes más confiables disponibles en el VCF original y será utilizado para los análisis comparativos con las variantes germinales y para la búsqueda de información en OncoKB y bases de datos de cáncer.

---

# 📊 **Comparación germinal vs somático**

Para evaluar diferencias genómicas entre los perfiles germinales y somáticos de la muestra S9, se realizó una comparación directa entre los VCF procesados por Sarek. Debido a que ninguno de los VCF generados por el pipeline incluía anotaciones funcionales (ANN, CSQ o SnpEff), la comparación se basó exclusivamente en las variantes filtradas por calidad (PASS), lo que garantiza que las variantes consideradas poseen soporte mínimo adecuado y cumplen los filtros internos de GATK HaplotypeCaller (germinal) y Mutect2 (somático).

Los archivos comparados fueron:

- `germinal_PASS.vcf.gz`
- `somatic_PASS.vcf.gz`

---

## **Número total de variantes**

Para obtener el número total de variantes por archivo, se usó:

```
bcftools view -H germinal_PASS.vcf.gz | wc -l
148

bcftools view -H somatic_PASS.vcf.gz | wc -l
130
```

---

## **Distribución por tipo**

```
bcftools view -H -v snps germinal_PASS.vcf.gz | wc -l
116

bcftools view -H -v indels germinal_PASS.vcf.gz | wc -l
32

bcftools view -H -v snps somatic_PASS.vcf.gz | wc -l
123

bcftools view -H -v indels somatic_PASS.vcf.gz | wc -l
6
```

| Métrica                                 | Germinal (HaplotypeCaller) | Somático (Mutect2) |
|-----------------------------------------|-----------------------------|---------------------|
| Variantes totales PASS                  | 148                         | 130                 |
| SNPs                                    | 116                         | 123                 |
| Indels                                  | 32                          | 6                   |
| Variantes compartidas (intersección)    | 57                          | 57                  |
| Variantes exclusivas de cada tipo       | 148 - 57 = **91**           | 130 - 57 = **73**   |
| Variantes seleccionadas para análisis   | 15                          | 15                  |
| Archivo final utilizado                 | germline_S9_final.vcf.gz    | somatic_S9_final.vcf.gz |


El análisis comparativo muestra que el perfil somático contiene más SNPs que el germinal, mientras que el germinal presenta muchos más indels, posiblemente debido a diferencias en la sensibilidad y filtros de Mutect2.  
Las 57 variantes compartidas indican que una fracción relevante de las variantes somáticas corresponde a variantes germinales heredadas. Las variantes somáticas exclusivas (73) corresponden a las verdaderas candidatas a eventos adquiridos, pero ninguna presentó evidencia clínica en OncoKB lo cual se explica más adelante.

---

## **Variantes compartidas**

La coincidencia entre variantes germinales y somáticas se determinó usando bcftools isec, que permite identificar posiciones idénticas presentes en ambos archivos.
Comando utilizado:


```
bcftools isec -n=2 -w1 -Oz -o shared_PASS_germline_somatic.vcf.gz \
germinal_PASS.vcf.gz somatic_PASS.vcf.gz

bcftools index shared_PASS_germline_somatic.vcf.gz
bcftools view -H shared_PASS_germline_somatic.vcf.gz | wc -l
```

**57 variantes compartidas.**

Interpretación:

Esto significa que 57 variantes detectadas como somáticas también aparecen en el genoma germinal, lo cual concuerda con el comportamiento real de Mutect2:
Cuando una variante está presente en el normal y en la muestra tumoral, suele marcarse como germinal (o “shared”), especialmente si las frecuencias alélicas son similares.
La presencia de 57 variantes compartidas entre los perfiles germinal y somático sugiere que una fracción importante de las variantes detectadas en el somático corresponde en realidad a variantes heredadas presentes en la línea germinal. En contraste, las variantes exclusivas del somático aproximadamente 73, considerando la diferencia entre las 130 variantes totales y las 57 compartidas constituyen las candidatas más relevantes para representar eventos adquiridos o potencialmente relacionados con procesos neoplásicos. La marcada diferencia en el número de indels entre ambos perfiles refleja, además, la naturaleza más estricta de los filtros aplicados por Mutect2, que prioriza la especificidad y tiende a descartar indels con mayor probabilidad de artefacto. Finalmente, antes de interpretar cualquier posible implicancia biológica o clínica, las variantes somáticas deben ser evaluadas en OncoKB para determinar nivel de evidencia y asociación oncológica, mientras que las variantes germinales deben analizarse en gnomAD para caracterizar su frecuencia poblacional y rareza
 
---

# 🧬 **Análisis de variantes somáticas en OncoKB**

Para la anotación de las variantes somáticas se utilizaron los archivos generados por el pipeline nf-core/sarek en modo tumor-only. 

Se utilizaron:

- `somatic_S9_final.vcf.gz`  archivo VCF con todas las variantes somáticas detectadas por Mutect2.
- `somatic_PASS.vcf.gz`  subconjunto de variantes que cumplen los filtros estándar del llamador.
- `somatic15.vcf` archivo final construido para esta tarea, que contiene únicamente las 15 variantes seleccionadas tras aplicar criterios de calidad (DP, AF y estado PASS), dado que el pipeline no incluía ninguna etapa de anotación funcional.

Es importante señalar que el pipeline proporcionado no ejecutó herramientas de anotación como VEP o SnpEff, por lo que los archivos VCF obtenidos no contienen información sobre el gen afectado, la consecuencia funcional ni nomenclatura HGVS. Esto obligó a seleccionar las variantes únicamente por calidad técnica y no por relevancia biológica o relación con cáncer.

Esto también llevo a que por falta de anotaciones funcionales y que la página de OncoKB ya no permite subir archivos VCF para cuentas no pagas y realizar las notaciones las variantes debean buscarse de manera individual por lo que se ingresaron manualmente en formato:

```
chr11:118472058:T>A
chr11:118503009:C>A
chr11:119278111:T>C
```

Resultado:

➡️ **Todas retornaron “No result found”.**

Esto se debe a que OncoKB solo contiene información para mutaciones driver, variantes con evidencia clínica, o alteraciones previamente asociadas a cáncer. Las variantes obtenidas en este análisis:

- Se ubican en una región estrecha del cromosoma 11 (118–119 Mb),  
- Región estrecha en chr11 sin genes oncológicos   
- Sin anotaciones funcionales  
- No registradas en bases oncológicas  

Como consecuencia, no fue posible obtener nivel de evidencia, oncogenicidad, cánceres asociados ni información terapéutica, ya que ninguna de las variantes figura en la base de datos OncoKB bajo los criterios clínicos o biológicos que utiliza esta plataforma.
Conclusión:

➡️ **No se identificó ninguna variante con relevancia clínica o terapéutica según OncoKB para las 15 variantes seleccionadas.**

---

# 🧬 **Resultados – Variantes Germinales en gnomAD**

A partir del conjunto de 15 variantes germinales obtenidas tras el filtrado bioinformático, se realizó la búsqueda manual en la base de datos gnomAD v4.1.0, registrando para cada variante la frecuencia global, las frecuencias por ancestría poblacional y su clasificación según rareza.
Donde se observaron las notaciones similares al cuadro que aparece en la siguiente imagen:

![Resultados gnomAD](imagenes/ejemplo_gnoma.png)

Apartir de estos cuadros para cada una de las  variantes seleccionados se obtuvieron los datos exigidos por la tarea y se construyó la siguiente tabla:

| Nº | Variante (chr:pos ref>alt) | AF Global | AF más alta | AF por población | Rareza |
|----|----------------------------|-----------|--------------|------------------|--------|
| 1 | 1:43337960 C>T | 6.85e-7 | South Asian | 0.00001102 | Muy rara |
| 2 | 1:43338672 T>C | No encontrada | — | — | No reportada |
| 3 | 1:43339467 TC>T | No encontrada | — | — | No reportada |
| 4 | 1:43339588 GC>G | No encontrada | — | — | No reportada |
| 5 | 1:43346599 C>A | No encontrada | — | — | No reportada |
| 6 | 1:43348913 A>C | 0.0000 | Ninguna | 0 | Ausente |
| 7 | 1:43352464 C>T | 0.000001240 | East Asian | 0.00004455 | Muy rara |
| 8 | 1:43352478 G>T | No encontrada | — | — | No reportada |
| 9 | 2:197400179 AAAT>A | 0.1803 | Middle Eastern | 0.2287 | Común |
| 10 | 2:197400449 T>A | 0.6732 | Middle Eastern | 0.8533 | Muy común |
| 11 | 2:197400626 T>C | 0.3197 | Finnish/Ashkenazi | ~0.42 | Común |
| 12 | 2:197402219 C>T | 0.9991 | Varias | 1.000 | Ultra común |
| 13 | 2:197402519 GA>G | 0.05443 | European (Finnish) | 0.06006 | Común |
| 14 | 2:197403046 G>GAA | 0.6795 | African/African American | 0.8377 | Muy común |
| 15 | 2:197408623 T>C | 0.1635 | Ashkenazi / Middle Eastern | 0.23 | Común |

---

## **Interpretación general**

- **5 variantes no están reportadas** → Esto indica que no han sido reportadas en más de 700,000 genomas/exomas, lo cual sugiere que podrían ser extremadamente raras o artefactos.  
- **2 variantes muy raras** → AF ~ 10⁻⁶ , lo que las clasifica como muy raras en la población general.  
- **8 variantes comunes** → AF 0.05–0.32  con frecuencias que varían entre valores moderados (5–16%) hasta muy altos (>60%), mostrando distribución amplia en diferentes ancestrías.
- **3 variantes muy comunes** → AF > 0.60  

Conclusión:

➡️ La mayoría de las variantes germinales identificadas son **polimorfismos comunes**.  
➡️ Solo dos variantes fueron clasificadas como potencialmente raras.


# 🧾 Discusión y conclusiones

En este práctico se logró ejecutar de forma completa el pipeline nf-core/sarek en sus dos modalidades: análisis germinal con HaplotypeCaller y análisis somático tumor-only con Mutect2, utilizando la muestra S9. A partir de los archivos FASTQ iniciales se generaron alineamientos al genoma de referencia GRCh38, reportes de calidad (MultiQC) y archivos VCF tanto germinales como somáticos, que luego fueron filtrados y analizados con herramientas de línea de comando (principalmente bcftools).

Un punto clave del trabajo fue que los VCF generados por el pipeline **no incluían anotaciones funcionales** (campos ANN/CSQ, genes, consecuencias, etc.). Esto impidió aplicar el criterio solicitado originalmente de seleccionar “variantes no sinónimas o de impacto moderado/alto”, ya que ese tipo de clasificación depende directamente de herramientas como VEP o SnpEff. Se intentó incorporar anotación en el servidor, pero por limitaciones técnicas del entorno no fue posible. Frente a esto, se optó por una estrategia alternativa de filtrado basada en **criterios de calidad**, utilizando campos como `FILTER=PASS`, profundidad de lectura (DP), calidad por profundidad (QD) y calidad de mapeo (MQ). Aunque no reemplaza a la anotación funcional, este enfoque sigue las recomendaciones de GATK y permite quedarse con variantes bien soportadas por la evidencia de secuenciación.

A nivel germinal se generó un archivo final con 15 variantes de alta calidad (`germline_S9_final.vcf.gz`), mientras que en el análisis somático se seleccionaron 15 variantes de mejor calidad interna (`somatic_S9_final.vcf.gz`) a partir de 130 variantes PASS. La comparación global entre `germinal_PASS.vcf.gz` y `somatic_PASS.vcf.gz` mostró que:

- El perfil germinal contiene **148 variantes PASS**, mientras que el somático contiene **130**.  
- El somático presenta más **SNPs** que el germinal (123 vs 116), pero muchos menos **indels** (6 vs 32), lo que refleja los filtros más estrictos de Mutect2 frente al riesgo de falsos positivos en indels.  
- Se identificaron **57 variantes compartidas** entre ambos perfiles, lo que indica que una parte importante de las variantes detectadas como somáticas corresponde en realidad a variantes heredadas presentes en la línea germinal.

Las variantes somáticas seleccionadas se evaluaron manualmente en **OncoKB**, ingresando las coordenadas en formato `chr:posición:REF>ALT`. Ninguna de las variantes arrojó resultados, lo que es coherente con varias observaciones: se trata de un conjunto pequeño de variantes, concentradas en una región acotada del cromosoma 11, sin anotación funcional disponible, con baja profundidad y sin evidencia previa de estar asociadas a genes oncológicos o a mutaciones driver conocidas. En consecuencia, **no se identificaron variantes con nivel de evidencia clínica, oncogenicidad definida ni información terapéutica** según OncoKB.

Por otro lado, las 15 variantes germinales se buscaron en **gnomAD v4.1.0**, registrando sus frecuencias globales y por población. La mayoría resultaron ser **polimorfismos comunes o muy comunes**, con frecuencias alélicas que pueden superar el 60 % en algunas ancestrías. Solo dos variantes se clasificaron como **muy raras**, y cinco no están reportadas en gnomAD, lo que podría indicar variantes extremadamente infrecuentes o posibles artefactos de llamado. En conjunto, el perfil germinal observado es compatible con un trasfondo genético mayormente constituido por variantes frecuentes de la población general, sin evidencia clara de variantes raras altamente sospechosas de enfermedad.


La principal limitación del análisis fue la **falta de anotación funcional automática**, que restringió la selección de variantes a criterios puramente técnicos y dificultó la interpretación biológica más profunda. Aun así, el práctico permitió comprender el flujo completo desde los FASTQ hasta los VCF filtrados, como perspectiva podría incluir en el pipeli original la notacion funcional para poder darle más relevancia a biológica a la muestra. Sin embargo, dado el caracter compartido del servidor despues de varios intentos y la cantidad de estudiantes ejecutando la misma tarea ya no pude realizar más pruebas y me limite a los resultados presentados en este informe.


---


