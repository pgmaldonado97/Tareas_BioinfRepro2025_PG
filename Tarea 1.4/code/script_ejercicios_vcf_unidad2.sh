#!/usr/bin/env bash
# Tarea 1.4 — Script para responder preguntas Unidad 2 seccion vcf
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p results data

# P1) ¿Cuántos individuos y variantes (SNPs) tiene el archivo?
zgrep -v "^#" /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf.gz | wc -l
zgrep -v "^#" /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf.gz | wc -l
zgrep -v "^#" /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf.gz | wc -w

# P2) ¿Cuántos sitios del archivo no tienen datos perdidos?
vcftools --gzvcf /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf.gz \
  --max-missing 1.0 \
  --out results/missing_site

# P3) Genera un archivo que contenga solo SNPs en una ventana de 2 Mb en cualquier cromosoma.
#     Nómbralo CLG_Chr<X>_<Start>-<End>Mb.vcf (usaste chr3: 30000–10000000)
vcftools --gzvcf /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf.gz \
  --chr 3 --from-bp 30000 --to-bp 10000000 \
  --recode --recode-INFO-all \
  --out results/CLG_Chr3_30000_10000000_Mb.vcf

# P4) Reporta cuántas variantes tiene el archivo generado.
grep -cv '^#' results/CLG_Chr3_30000_10000000_Mb.vcf.recode.vcf

# P5) Reporta la cobertura promedio para todos los individuos del set de datos.
cd results
vcftools --vcf CLG_Chr3_30000_10000000_Mb.vcf --depth --out depth
head -n 20 depth.idepth || true

# P6) Calcula la frecuencia de cada alelo para todos los individuos dentro del archivo y guarda el resultado.
vcftools --vcf CLG_Chr3_30000_10000000_Mb.vcf --freq --out allele_freq

# P7) Filtra el archivo de frecuencias para solo incluir variantes bialélicas.
awk '$3 == 2' allele_freq.frq > allele_freq_biallelic.frq
# less allele_freq_biallelic.frq

# P8) Llama a un script en R que lee el archivo de frecuencias bialélicas y guarda un histograma con el espectro de MAF.
cd ..
cat > code/script_R_ejercicio8.R <<'RSCRIPT'
infile  <- "results/allele_freq_biallelic.frq"
out_tsv <- "results/maf_values.tsv"
out_png <- "results/maf_histogram.png"

df <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
f1 <- as.numeric(sub(".*:", "", df[[5]]))
f2 <- as.numeric(sub(".*:", "", df[[6]]))
maf <- pmin(f1, f2, na.rm = TRUE)

out <- data.frame(MAF = maf)
write.table(out, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

png(out_png, width = 1200, height = 800, res = 150)
hist(maf, breaks = 50, main = "MAF spectrum (biallelic)", xlab = "Minor allele frequency")
dev.off()

cat("Script completado. Archivos generados:\n- ", out_tsv, "\n- ", out_png, "\n", sep = "")
RSCRIPT
Rscript code/script_R_ejercicio8.R

# P9) ¿Cuántos sitios tienen una frecuencia del alelo menor <0.05?
cd results
awk -F'\t' 'NR>1 && $1 < 0.05 {n++} END{print n+0}' maf_values.tsv

# P10) Calcula la heterocigosidad de cada individuo.
vcftools --vcf CLG_Chr3_30000_10000000_Mb_biallelic.vcf --het --out heterozygosity
ls
head heterozygosity.het

# P11) Calcula la diversidad nucleotídica por sitio.
vcftools --vcf CLG_Chr3_30000_10000000_Mb.vcf --site-pi --out nucleotide_diversity
ls
head -n 10 nucleotide_diversity.sites.pi

# P12) Filtra los sitios que tengan una frecuencia del alelo menor <0.05.
vcftools --vcf CLG_Chr3_30000_10000000_Mb_biallelic.vcf --max-maf 0.05 --recode --out allele_biallelic_maf05
ls
grep -vc '^#' allele_biallelic_maf05.recode.vcf

# P13) Convierte el archivo a formato PLINK.
 /opt/plink/plink-1.90/plink \
  --vcf allele_biallelic_maf05.recode.vcf \
  --double-id \
  --allow-extra-chr \
  --snps-only just-acgt \
  --make-bed \
  --out allele_biallelic_maf05
ls

