# Tarea 1.4 — Informe
#(Ejercicio 4)

**Directorio de trabajo:** `Tarea 1.4/`  

## Estructura
~~~
Tarea 1.4/
├── code/      # scripts (R y bash)
├── data/      # entradas: PLINK + metadatos
└── results/   # salidas intermedias y finales
~~~

## Software
- R 4.4.3 (`/opt/R/R-4.4.3/bin/Rscript`)
- PLINK 1.9 (`/opt/plink/plink-1.90/plink`)

---

## Objetivo del Ejercicio 4

Usar `data/chilean_all48_hg19_popinfo.csv` y el comando `plink --update-ids` para **cambiar el FID** de las muestras del set PLINK `data/chilean_all48_hg19.{bed,bim,fam}`, de manera que el **FID** pase a ser **`Categ.Altitud`** tomado desde el metadato de maíz `data/maizteocintle_SNP50k_meta_extended.txt`.

> Para que esto sea posible debe existir una **llave en común** (columna con valores coincidentes) entre las 48 muestras de `popinfo` y las filas relevantes del meta de maíces; sin una llave, no se puede construir la tabla requerida por `--update-ids`:
> ~~~
> oldFID  oldIID  newFID  newIID
> ~~~

---

## Entradas utilizadas
- `data/chilean_all48_hg19.{bed,bim,fam}`
- `data/chilean_all48_hg19_popinfo.csv`
- `data/maizteocintle_SNP50k_meta_extended.txt`

---

## Pasos realizados

### 1) Extraer IDs del `.fam`

**Script:** `code/step4a_extract_fam_ids.R`  
**Comando:**
~~~bash
/opt/R/R-4.4.3/bin/Rscript code/step4a_extract_fam_ids.R
~~~

**Salida principal**
- `results/fam_ids.tsv` — columnas: `oldFID`, `oldIID`, `IID_STD` (ID estandarizado).

> *Nota:* Este paso se usó para comparar correctamente los IIDs del `.fam` con campos de `popinfo` y asegurar que estamos trabajando sobre las mismas 48 muestras.

---

### 2) Escaneo de llaves entre `popinfo` y el meta

Para **no asumir** cuál es la llave, se escanearon **todas** las columnas de `popinfo` contra **todas** las columnas del archivo `maizteocintle_SNP50k_meta_extended.txt`, contando coincidencias tras una normalización (mayúsculas, `trim`, y eliminación de acentos/ñ).

**Script:** `code/step4z_scan_keys.R`  
**Comando:**
~~~bash
/opt/R/R-4.4.3/bin/Rscript code/step4z_scan_keys.R
~~~

**Salidas**
- `results/key_scan_fam_pop.tsv` — confirma que **`IndID`** de `popinfo` **matchea 48/48** con los IIDs del `.fam`.
- `results/key_scan_pop_meta.tsv` — **todas** las combinaciones `pop_col × meta_col` tienen **`matches = 0`**.

**Interpretación**
- Aunque `IndID` enlaza perfecto `fam ↔ popinfo`, **no existe llave** con el meta de maíces.  
- Columnas candidatas como `Num_Colecta`, `NSiembra`, `OrderColecta`, `Categ.Altitud`, etc., **no presentan intersección de valores** con estas 48 muestras.

---

### 3) Consecuencia sobre `plink --update-ids`

Dado que **no hay llave en común**, **no se puede** construir la tabla de 4 columnas `oldFID oldIID newFID newIID` con `newFID = Categ.Altitud` del meta.  
Por lo tanto, **no se ejecutó** el siguiente comando:

~~~bash
/opt/plink/plink-1.90/plink \
  --bfile data/chilean_all48_hg19 \
  --update-ids results/update_ids.txt \
  --make-bed \
  --out results/chilean_all48_hg19_altfid
~~~

---

## Evidencia

1) **`results/key_scan_fam_pop.tsv`**  
Confirma que `IndID` (en `popinfo`) coincide **48/48** con los IIDs del `.fam`.  
➡️ `popinfo` corresponde a las mismas 48 muestras del `.fam`.

2) **`results/key_scan_pop_meta.tsv`**  
Todas las filas muestran **`matches = 0`** para cualquier par `pop_col × meta_col`.  
➡️ **No existe llave** entre `popinfo` (48 chilenos) y `maizteocintle_SNP50k_meta_extended.txt`.

---

## Conclusión del Ejercicio 4

Con el escaneo exhaustivo de llaves, se concluye que **no existe una llave común** entre las 48 muestras chilenas y el metadato de maíces. En consecuencia, **no es posible** renombrar el FID a `Categ.Altitud` usando `plink --update-ids` para este conjunto de datos.  
Se documenta el intento, los scripts utilizados y la evidencia generada. Se continúa con los **Ejercicios 5–7**, que solo requieren `fam` y `popinfo`.

> **Plan B (no aplicado para respetar el enunciado):** podría derivarse una categoría de altitud **directamente desde `popinfo`** (si existe una columna numérica o categórica equivalente) y usarla en `--update-ids`. En este informe se **opta por no hacerlo** porque la consigna exige explícitamente **`Categ.Altitud` del meta** de maíces.

---

## Reproducibilidad rápida

Desde `Tarea 1.4/`:

~~~bash
# 1) Extraer IDs del fam
/opt/R/R-4.4.3/bin/Rscript code/step4a_extract_fam_ids.R

# 2) Escanear llaves popinfo × meta
/opt/R/R-4.4.3/bin/Rscript code/step4z_scan_keys.R
~~~

**Revisar:**
- `results/fam_ids.tsv`
- `results/key_scan_fam_pop.tsv`
- `results/key_scan_pop_meta.tsv`  ← aquí se observa `matches = 0` para todas las parejas.

---

## Scripts y salidas (Ejercicio 4)

**Scripts**
- `code/step4a_extract_fam_ids.R`
- `code/step4z_scan_keys.R`

**Salidas**
- `results/fam_ids.tsv`
- `results/key_scan_fam_pop.tsv`
- `results/key_scan_pop_meta.tsv`

**Entradas**
- `data/chilean_all48_hg19.{bed,bim,fam}`
- `data/chilean_all48_hg19_popinfo.csv`
- `data/maizteocintle_SNP50k_meta_extended.txt`

