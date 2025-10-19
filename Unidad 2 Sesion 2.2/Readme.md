# Unidad 2 Sesion 2.2

Esta tarea corresponde a la Sesion 2 Unidad 2 titulado **Análisis genético de poblaciones**

Para hacer segumiento a los comando puede usarse el repositorio del curso:

`Unidad2/Sesion2/Tutorial_PopGeno.md`

Donde se detalla todos los comandos y explicaciones necesarias del análisis.

---

## Paso 1: Investigar la cantidad de genotipos perdidos los datos de 1000G

### ¿Cómo se llaman los archivos que contienen las tasas de datos perdidos por SNP y por muestra?

PLINK siempre genera:

- `nombre.imiss` → por muestra → `chilean_all48_hg19.imiss`
- `nombre.lmiss` → por SNP → `chilean_all48_hg19.lmiss`

### ¿Cuántas variantes se eliminaron por tener una tasa de datos perdidos mayor a 0.2?

4680 es el número de SNPs removidos

```bash
grep -i -E 'removed.*missing genotype data' chilean_all48_hg19_2.log | awk '{print $1}'
```

### ¿Cuántos individuos tenían una tasa de datos perdidos mayor a 0.02?

```bash
awk 'NR>1 && $6>0.02 {c++} END{print c+0}' plink.imiss
```

Dos muestras removidas

### Basados en los histogramas y en sus cálculos, ¿qué valores umbrales de datos perdidos para muestras y SNPs sugeriría?

Con base en tus histogramas:

**Umbral para SNPs (–geno):**

Lo más recomendado es 0.02 (2%), porque la mayoría de SNPs están muy bien y solo unos pocos tienen valores más altos.

Si quieres ser más permisiva y no perder demasiados SNPs, se podría usar 0.05 (5%), pero 0.02 es un buen estándar.

**Umbral para individuos (–mind):**

El histograma muestra que la mayoría de las muestras tienen <0.02, pero unas pocas tienen 0.10–0.15.

Con un corte de 0.02, eliminarías justamente a esas muestras malas y te quedas con individuos de alta calidad.

Esto coincide con lo que comentaste que pasó en clases: al bajar el umbral, sí se eliminaban ciertos individuos.

---

## Paso 2: Revisar discrepancias entre el sexo declarado en la tabla de fenotipos (archivo ped/bed) y el inferido desde los genotipos.

### ¿Cuántos individuos fueron eliminados por discrepancia de sexo?

```bash
wc -l sex_discrepancy.txt
```

3 fueron eliminados

### ¿Qué riesgo(s) se corre(n) si no se eliminaran?

Si no eliminas los individuos con discrepancia de sexo:

- **Error de anotación fenotípica** → el sexo declarado no coincide con el sexo genotípico real.  
- **Sesgo en análisis posteriores** → como estudios de asociación, imputación, o análisis de cromosoma X, ya que esos individuos no cumplen los supuestos.  
- **Riesgo de intercambios de muestra** → puede indicar que hubo un *swap* (ej. muestra de varón asignada a una mujer).  
- **Impacto en control de calidad** → se mantienen datos posiblemente contaminados o mal rotulados que distorsionan la interpretación.

---

## Paso 3: Filtrado de SNPs

### ¿Cuál es el nombre del primer conjunto de datos que solo contiene SNPs en autosomas?

El nombre es `chilean_all48_hg19_7`.

### ¿Cuántos SNPs se encontraban en cromosomas sexuales?

```bash
expr $(wc -l < chilean_all48_hg19_6.bim) - $(wc -l < chilean_all48_hg19_7.bim)
```

16,702 SNP

### ¿Cómo calcularía el número de cromosomas que porta cada uno de los alelos para cada SNP?

En la práctica, PLINK lo hace automáticamente con la opción `--freq counts`, que entrega el número exacto de cromosomas portadores de cada alelo (columnas C1 y C2 en el archivo `allele_counts.frq.counts`).

```bash
plink --bfile chilean_all48_hg19_7 --freq counts --out allele_counts
```

---

## Paso 4: Borrar SNPs por filtro de HWE

### ¿Cuál es el nombre del archivo con los resultados de la prueba de HWE?

`plink.hwe`

### ¿Basándose en la distribución de los valores de p, le parece el umbral usado razonable o propondría otro valor?

El umbral de **1e-6** es razonable, porque filtra solamente a los SNPs con desviaciones extremas (los que aparecen en el segundo histograma).

Mantener este umbral ayuda a eliminar variantes probablemente espurias sin perder la gran mayoría de SNPs que cumplen bien con HWE.

Podría usarse un umbral un poco menos estricto (**1e-5**) si se quisiera ser más permisivo, pero **1e-6** es un estándar recomendado en control de calidad de GWAS y se ajusta bien a los datos que muestran tus histogramas.

---

## Paso 5: Eliminar parentescos desconocidos

### ¿Cuántos SNPs en aparente equilibrio de ligamiento se encontraron?

En primer lugar, eliminamos los SNPs que están muy correlacionados entre sí (es decir, en desequilibrio de ligamiento) y conserva solo los independientes.

```bash
plink --bfile chilean_all48_hg19_9 --exclude $T/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
```

Los SNPs que quedaron en “equilibrio de ligamiento” se guardan en el archivo: `indepSNP.prune.in`

Y los que se eliminaron (porque estaban en fuerte LD) se guardan en: `indepSNP.prune.out`

Y finalmente contamos de aquí el número:

```bash
wc -l indepSNP.prune.in
```

De eso se obtiene que son **103,214** los SNPs independientes que quedaron.

### ¿Cuántos SNPs se eliminaron por estar en regiones de inversiones conocidas?

El archivo `inversion.txt` tiene las coordenadas de regiones de inversión conocidas (zonas de cromosomas con fuerte LD por razones estructurales).  
Entonces al contar cuántos SNPs fueron excluidos de esas regiones podemos obtener la respuesta.  
Eso lo logramos con el siguiente comando:

```bash
wc -l indepSNP.prune.out
```

El número obtenido fue de **346,968**.

### ¿Cuántos individuos quedaron luego del filtro de parentesco?

Primero utilizamos el comando:

```bash
plink --bfile chilean_all48_hg19_9 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
```

Y luego contamos cuántos individuos quedaron con el siguiente comando:

```bash
wc -l chilean_all48_hg19_10.fam
```

El número obtenido fue de **42 individuos**.

### ¿Cuál fue el mayor coeficiente de parentesco efectivamente aceptado?

El **PI_HAT** mide la proporción de genoma compartido entre dos individuos.

Umbral: se eliminó si PI_HAT ≥ 0.2.

Así que el mayor coeficiente aceptado será el valor máximo justo por debajo de 0.2.

Tras excluir pares con PI_HAT ≥ 0.2, el mayor parentesco aceptado fue **PI_HAT = 0.1996**, obtenido usando los siguientes comandos:

```bash
plink --bfile chilean_all48_hg19_9 --extract indepSNP.prune.in --genome --out pihat_all
awk 'NR>1 && $10<0.2 {if($10>m)m=$10} END{print m+0}' pihat_all.genome
```

---

## Paso 6: Graficar resultados de MDS

En R, genere gráficos similares para las combinaciones **Component 2 vs 3** y **3 vs 4**. ¿Qué puede concluir de estos gráficos?

Primero generé un script en R para graficar las combinaciones solicitadas de componentes principales (2 vs 3 y 3 vs 4) a partir de los resultados del MDS.  
El script carga los datos de las coordenadas (`MDS_merge2.mds`) y la información de etnicidad de cada individuo (`ethnicityfile.txt`), los combina y finalmente genera dos gráficos en formato PNG.

### Script en R

```r
mds <- read.table("MDS_merge2.mds", header = TRUE)
eth <- read.table("ethnicityfile.txt", header = TRUE)
data <- merge(mds, eth, by.x = "IID", by.y = "IID")

png("MDS_Comp2_vs_Comp3.png", width=800, height=600)
plot(data$C2, data$C3, col=as.factor(data$ethnicity), pch=19,
     xlab="MDS Component 2", ylab="MDS Component 3",
     main="MDS: Component 2 vs Component 3")
legend("topright", legend=levels(as.factor(data$ethnicity)),
       col=1:length(unique(data$ethnicity)), pch=19)
dev.off()

png("MDS_Comp3_vs_Comp4.png", width=800, height=600)
plot(data$C3, data$C4, col=as.factor(data$ethnicity), pch=19,
     xlab="MDS Component 3", ylab="MDS Component 4",
     main="MDS: Component 3 vs Component 4")
legend("topright", legend=levels(as.factor(data$ethnicity)),
       col=1:length(unique(data$ethnicity)), pch=19)
dev.off()
```

Las imágenes generadas se pueden observar en las figuras:

- `MDS_Comp2_vs_Comp3.png`
- `MDS_Comp3_vs_Comp4.png`

---

## Paso 7: Realizar un análisis de ancestría

### ¿Cuántos SNPs quedaron luego del filtro?

Para determinarlo no fue necesario ejecutar **ADMIXTURE** completo, ya que este análisis es muy demandante en términos computacionales.  
En su lugar, se utilizó el archivo de salida generado con el comando de *LD pruning* (`MDS_merge_r2_lt_0.2.bim`) y se contó el número de filas, cada una correspondiente a un SNP:

```bash
wc -l MDS_merge_r2_lt_0.2.bim
```

El resultado obtenido fue de **70,534 SNPs** retenidos después del filtro.

### ADMIXTURE asume que los individuos no están emparentados. Sin embargo, no realizamos ningún filtro. ¿Por qué?

En pasos previos del pipeline ya habíamos eliminado a los individuos con parentescos no deseados (pihat ≥ 0.2), lo que asegura que la muestra final no incluye individuos estrechamente relacionados.  
Por esta razón no fue necesario aplicar un nuevo filtro en esta etapa, ya que el dataset que se usará en ADMIXTURE ya cumple con esta condición.


