# step4a_extract_fam_ids.R
# Lee el .fam y escribe un TSV con: oldFID, oldIID, IID_STD (Del ID estandarizado)
options(stringsAsFactors = FALSE)

fam_path <- "../data/chilean_all48_hg19.fam"
out_path <- "../results/fam_ids.tsv"

std <- function(x) trimws(toupper(as.character(x)))  # mayúsculas y sin espacios extremos

# El .fam tiene 6 columnas sin cabecera
fam <- read.table(fam_path, header = FALSE)
names(fam) <- c("FID","IID","PID","MID","SEX","PHENO")

res <- data.frame(
  oldFID  = fam$FID,
  oldIID  = fam$IID,
  IID_STD = std(fam$IID),
  stringsAsFactors = FALSE
)

dir.create("../results", showWarnings = FALSE, recursive = TRUE)
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Hecho.\n")
cat("Archivo generado:", out_path, "\n")
cat("Filas:", nrow(res), "\n")

