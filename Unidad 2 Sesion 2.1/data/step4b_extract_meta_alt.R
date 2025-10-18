# step4b_extract_meta_alt.R (versión con puente robusto)
options(stringsAsFactors = FALSE)

fam_ids_path <- "../results/fam_ids.tsv"
pop_path     <- "../data/chilean_all48_hg19_popinfo.csv"
meta_path    <- "../data/maizteocintle_SNP50k_meta_extended.txt"
out_path     <- "../results/meta_alt.tsv"

std <- function(x) {
  z <- as.character(x)
  z <- iconv(z, from = "", to = "ASCII//TRANSLIT", sub = "")
  toupper(trimws(z))
}

# 1) fam_ids
fam_ids <- read.delim(fam_ids_path, check.names = FALSE)
fam_ids$IID_STD <- std(fam_ids$IID_STD)

# 2) popinfo
pop <- tryCatch(read.csv(pop_path, check.names = FALSE),
                error=function(e) read.delim(pop_path, sep=",", header=TRUE, check.names=FALSE))
for(nm in names(pop)) if(!is.character(pop[[nm]])) pop[[nm]] <- as.character(pop[[nm]])

# 3) meta (manejo de codificación)
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

# 4) ¿qué columna de popinfo coincide con los IID del fam?
fam_set <- unique(fam_ids$IID_STD)
pop_match_counts <- sapply(pop, function(col) sum(std(col) %in% fam_set))
pop_col_fam <- names(pop_match_counts)[which.max(pop_match_counts)]
cat("Columna en popinfo que coincide con IID del fam:", pop_col_fam, "\n")

pop$POP_FAM_STD <- std(pop[[pop_col_fam]])

# 5) localizar Categ.Altitud en meta
alt_cols <- grep("^Categ(\\.|_|\\s)?Altitud$", names(meta), ignore.case = TRUE, value = TRUE)
if (!length(alt_cols)) stop("No encontré columna 'Categ.Altitud' en el meta.")
alt_col <- alt_cols[1]
meta$ALT <- as.character(meta[[alt_col]])

# 6) elegir puente popinfo -> meta
#    probamos candidatos habituales; si no hay, buscamos el mejor por scoring pero exigimos >0 matches
candidatos_pop  <- c("Num_Colecta","NSiembra","OrdenColecta","OrderColecta")
candidatos_meta <- c("Num_Colecta","NSiembra","OrdenColecta","OrderColecta")

# primero intenta pares con mismo nombre en ambos:
puente_ok <- FALSE
for(nm in candidatos_pop){
  if(nm %in% names(pop) && nm %in% names(meta)){
    pop_vec  <- unique(std(pop[[nm]]))
    meta_vec <- unique(std(meta[[nm]]))
    sc <- sum(pop_vec %in% meta_vec)
    if(sc > 0){
      pop_col_bridge <- nm
      meta_col_id    <- nm
      puente_ok <- TRUE
      cat("Puente (candidato directo) elegido:", nm, "->", nm, "(coincidencias:", sc, ")\n")
      break
    }
  }
}

# si lo anterior no funcionó, buscar el mejor par por scoring y exigir >0
if(!puente_ok){
  best <- list(score=-1, pop_col=NA, meta_col=NA)
  for(pc in names(pop)){
    pj <- unique(std(pop[[pc]]))
    if(all(!nzchar(pj))) next
    scores <- sapply(names(meta), function(mc) sum(pj %in% unique(std(meta[[mc]]))))
    sc <- max(scores, na.rm=TRUE)
    if(length(sc)==0 || is.infinite(sc)) sc <- 0
    if(sc > best$score){
      best$score <- sc
      best$pop_col <- pc
      best$meta_col <- names(meta)[which.max(scores)]
    }
  }
  if(best$score > 0){
    pop_col_bridge <- best$pop_col
    meta_col_id    <- best$meta_col
    puente_ok <- TRUE
    cat("Puente (auto) elegido:", pop_col_bridge, "->", meta_col_id, "(coincidencias:", best$score, ")\n")
  }
}

if(!puente_ok){
  stop("No se encontró un par de columnas puente con coincidencias > 0.\n",
       "Sugerencia: revisa nombres en popinfo y meta y dime cuál usar (p.ej. Num_Colecta).")
}

# 7) estandarizar columnas puente y hacer merges
pop$POP_BRIDGE_STD <- std(pop[[pop_col_bridge]])
meta$META_ID_STD   <- std(meta[[meta_col_id]])

m1 <- merge(fam_ids[, "IID_STD", drop=FALSE],
            pop[, c("POP_FAM_STD","POP_BRIDGE_STD"), drop=FALSE],
            by.x="IID_STD", by.y="POP_FAM_STD", all.x=TRUE)

m2 <- merge(m1, meta[, c("META_ID_STD","ALT"), drop=FALSE],
            by.x="POP_BRIDGE_STD", by.y="META_ID_STD", all.x=TRUE)

res <- data.frame(IID_STD = m2$IID_STD, ALT = m2$ALT, stringsAsFactors = FALSE)
keep <- !is.na(res$ALT) & nzchar(res$ALT) & res$ALT != "NA"
res  <- res[keep, , drop = FALSE]

dir.create("../results", showWarnings = FALSE, recursive = TRUE)
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Hecho.\nArchivo generado:", out_path, "\nFilas:", nrow(res), "\n")

