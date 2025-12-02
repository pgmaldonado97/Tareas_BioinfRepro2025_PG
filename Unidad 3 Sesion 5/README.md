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

El análisis bioinformático se realizó utilizando el pipeline **nf-core/sarek**, ejecutado mediante Nextflow y Singularity en un entorno Linux del servidor institucional.  
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

---

### **3. Ejecución del análisis germinal**

Se utilizó el script **sarek_germinal.sh**, basado en el modelo entregado y adaptado para la muestra S9.

Comando ejecutado:

```
bash sarek_germinal.sh S9_R1.fastq.gz S9_R2.fastq.gz ../results
```

Esto generó:

- VCF germinales  
- Reportes de calidad  
- Archivos BAM procesados  

---

### **4. Ejecución del análisis somático**

Se utilizó el script **sarek_somatic.sh**, que ejecuta Mutect2 en modo tumor-only.

Comando ejecutado:

```
bash sarek_somatic.sh S9_R1.fastq.gz S9_R2.fastq.gz ../results
```

El pipeline produjo los VCF somáticos y reportes de calidad.

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

---

# 🔬 **ANÁLISIS GERMINAL**

## **Filtrado de variantes germinales**

El archivo inicial fue:

```
variant_calling/haplotypecaller/S9/S9.haplotypecaller.filtered.vcf.gz
```

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

---

### **2. Estrategia alternativa de filtrado**

Dado lo anterior, se aplicó un filtrado basado en calidad:

- `FILTER = PASS`  
- QD alto  
- MQ > 40  
- DP > 10  

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

```
germline_S9_final.vcf.gz
```

Contiene 15 variantes germinales de alta calidad.

---

# 🔬 **Filtrado de Variantes Somáticas (Mutect2)**

Archivo inicial:

```
variant_calling/mutect2/S9/S9.mutect2.filtered.vcf.gz
```

---

### **1. Verificación de anotaciones**

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

---

# 📊 **Comparación germinal vs somático**

Los archivos comparados fueron:

- `germinal_PASS.vcf.gz`
- `somatic_PASS.vcf.gz`

---

## **Número total de variantes**

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

---

## **Variantes compartidas**

```
bcftools isec -n=2 -w1 -Oz -o shared_PASS_germline_somatic.vcf.gz \
germinal_PASS.vcf.gz somatic_PASS.vcf.gz

bcftools index shared_PASS_germline_somatic.vcf.gz
bcftools view -H shared_PASS_germline_somatic.vcf.gz | wc -l
```

**57 variantes compartidas.**

Interpretación:

- mutaciones compartidas → variantes germinales también vistas por Mutect2  
- somático exclusivo (~73 variantes) → candidatas somáticas reales  

---

# 🧬 **Análisis de variantes somáticas en OncoKB**

Se utilizaron:

- `somatic_S9_final.vcf.gz`  
- `somatic_PASS.vcf.gz`  
- `somatic15.vcf`

Dado que los VCF **no poseen anotación funcional**, las variantes se ingresaron manualmente en formato:

```
chr11:118472058:T>A
chr11:118503009:C>A
chr11:119278111:T>C
```

Resultado:

➡️ **Todas retornaron “No result found”.**

Razones:

- Ninguna corresponde a variantes driver  
- Región estrecha en chr11 sin genes oncológicos  
- Baja profundidad  
- Sin anotaciones funcionales  
- No registradas en bases oncológicas  

Conclusión:

➡️ **No se identificó ninguna variante con relevancia clínica o terapéutica según OncoKB.**

---

# 🧬 **Resultados – Variantes Germinales en gnomAD**

Se evaluaron las 15 variantes germinales en gnomAD v4.1.0.

Los resultados se registraron en esta tabla:

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

- **5 variantes no están reportadas** → extremadamente raras o artefactos.  
- **2 variantes muy raras** → AF ~ 10⁻⁶  
- **8 variantes comunes** → AF 0.05–0.32  
- **3 variantes muy comunes** → AF > 0.60  

Conclusión:

➡️ La mayoría de las variantes germinales identificadas son **polimorfismos comunes**.  
➡️ Solo dos variantes fueron clasificadas como potencialmente raras.

---


