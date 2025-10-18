# step4z_scan_keys.R — Diagnóstico de llaves popinfo ↔ meta y fam ↔ popinfo
options(stringsAsFactors = FALSE)

fam_ids_path <- "../results/fam_ids.tsv"
pop_path     <- "../data/chilean_all48_hg19_popinfo.csv"
meta_path    <- "../data/maizteocintle_SNP50k_meta_extended.txt"
out_pop_meta <- "../results/key_scan_pop_meta.tsv"
out_fam_pop  <- "../results/key_scan_fam_pop.tsv"

std <- function(x){
  z <- as.character(x)
  z <- iconv(z, from = "", to = "ASCII//TRANSLIT", sub = "") # quita acentos/ñ
  toupper(trimws(z))
}

# --- Carga ---
fam <- read.delim(fam_ids_path, check.names = FALSE)

pop <- tryCatch(read.csv(pop_path, check.names = FALSE),
                error=function(e) read.delim(pop_path, sep=",", header=TRUE, check.names=FALSE))
for(nm in names(pop)) if(!is.character(pop[[nm]])) pop[[nm]] <- as.character(pop[[nm]])

meta <- tryCatch(
  read.delim(meta_path, check.names = FALSE, fileEncoding = "UTF-8"),
  error = function(e) {
    tryCatch(
      read.delim(meta_path, check.names = FALSE, fileEncoding = "latin1"),
      error = function(e2) read.table(meta_path, header = TRUE, sep = "", check.names = FALSE)
    )
  }
)
for(nm in names(meta)) if(!is.character(meta[[nm]])) meta[[nm]] <- as.character(meta[[nm]])

# --- 1) ¿Qué columna de popinfo coincide mejor con los IID del fam? ---
fam_iids <- unique(std(fam$IID_STD))
fam_vs_pop <- data.frame(pop_col=character(), matches=integer(), stringsAsFactors=FALSE)

for(pc in names(pop)){
  pv <- unique(std(pop[[pc]]))
  if(all(!nzchar(pv))) next
  m <- sum(pv %in% fam_iids)
  fam_vs_pop <- rbind(fam_vs_pop, data.frame(pop_col=pc, matches=m))
}
fam_vs_pop <- fam_vs_pop[order(-fam_vs_pop$matches), ]
dir.create("../results", showWarnings = FALSE, recursive = TRUE)
write.table(fam_vs_pop, out_fam_pop, sep="\t", quote=FALSE, row.names=FALSE)
cat("Top columnas popinfo que matchean con IID del fam (ver", out_fam_pop, "):\n")
print(head(fam_vs_pop, 10))

# --- 2) Escaneo popinfo × meta ---
res <- data.frame(pop_col=character(), meta_col=character(),
                  matches=integer(), pop_uniq=integer(), meta_uniq=integer(),
                  example=character(), stringsAsFactors=FALSE)

# Precompute std uniques
pop_std <- lapply(pop, function(col) unique(std(col)))
meta_std <- lapply(meta, function(col) unique(std(col)))

for(pc in names(pop_std)){
  pv <- pop_std[[pc]]
  if(all(!nzchar(pv))) next
  for(mc in names(meta_std)){
    mv <- meta_std[[mc]]
    if(all(!nzchar(mv))) next
    inter <- intersect(pv, mv)
    m <- length(inter)
    if(m>0){
      ex <- paste(head(inter, 5), collapse=",")
    } else {
      ex <- ""
    }
    res <- rbind(res, data.frame(pop_col=pc, meta_col=mc,
                                 matches=m, pop_uniq=length(pv), meta_uniq=length(mv),
                                 example=ex, stringsAsFactors=FALSE))
  }
}

if(nrow(res)==0){
  cat("\nNo se encontraron coincidencias en NINGUNA pareja popinfo×meta.\n",
      "Conclusión: no hay llave común. Documenta esto y pasa a los ítems 5–7.\n")
} else {
  res <- res[order(-res$matches, res$pop_col, res$meta_col), ]
  write.table(res, out_pop_meta, sep="\t", quote=FALSE, row.names=FALSE)
  cat("\nTop pares popinfo×meta por coincidencias (ver", out_pop_meta, "):\n")
  print(head(res, 20))
}

