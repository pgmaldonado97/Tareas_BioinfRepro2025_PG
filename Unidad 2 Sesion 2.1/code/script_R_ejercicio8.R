# -------------------------------
# script_R_ejercicio8.R
# Este script lee el archivo de frecuencias bialélicas (.frq),
# calcula el MAF (minor allele frequency) y genera un histograma.
# -------------------------------

# 1. Definir rutas de entrada y salida
input_file  <- "results/allele_freq_biallelic.frq"    # archivo que cree en el ejercicio 7
output_tsv  <- "results/maf_values.tsv"               # archivo donde guardaremos los valores MAF
output_plot <- "results/maf_histogram.png"            # histograma en formato PNG

# 2. Leer el archivo .frq sin cabecera
frq <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)

# 3. Extraer las frecuencias de los dos alelos (últimas dos columnas)
get_freq <- function(x) { as.numeric(sub(".*:", "", x)) }
f1 <- get_freq(frq[[ncol(frq) - 1]])
f2 <- get_freq(frq[[ncol(frq)]])

# 4. Calcular el MAF como el mínimo entre las dos frecuencias
maf <- pmin(f1, f2)

# 5. Guardar los valores MAF en un archivo
write.table(maf, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = "MAF")

# 6. Crear un histograma del espectro de frecuencias alélicas
png(output_plot, width = 1000, height = 800, res = 120)
hist(maf,
     breaks = seq(0, 0.5, by = 0.01),
     main = "Espectro de frecuencias alélicas (MAF)",
     xlab = "MAF",
     ylab = "Número de SNPs")
dev.off()

# 7. Mensaje final
cat("Script completado. Archivos generados:\n")
cat("-", output_tsv, "\n")
cat("-", output_plot, "\n")

