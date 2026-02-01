# NA

```` markdown
abyssEdge
=========

Paquete R para detectar y analizar diferencias en expresión/transcripción y calcular métricas de diversidad (incluye funciones para Tsallis, pruebas estadísticas y visualización).

Nota sobre el origen
--------------------

abyssEdge hereda y adapta gran parte del código original de la librería SplicingFactory (uso y crédito al proyecto original), integrando además nuevas funciones, correcciones y utilidades específicas para el análisis de diversidad transcriptómica mediante la entropía de Tsallis.

Teoría de Tsallis
-------------------------

Tsallis es una generalización de la entropía de Shannon definida por
S_q = (1 - sum p^q) / (q - 1), donde p son las proporciones de isoformas.
Para q → 1 el resultado converge a la entropía de Shannon. En abyssEdge
se utiliza la entropía de Tsallis y sus correspondientes números de Hill (D_q)
para cuantificar la diversidad de isoformas por gen y compararla entre grupos
de muestras. El parámetro `q` permite sintonizar la sensibilidad a especies
raras o abundantes (valores q < 1 favorecen rarezas, q > 1 favorecen
abundantes).

Funcionalidades integradas
-------------------------

- Cálculo de entropía y diversidad por gen (Tsallis):
  - `calculate_tsallis_entropy`: calcula S_q y/o D_q para un vector de expresiones (soporta normalización, múltiples valores de `q` y el límite q→1). Devuelve vectores numéricos o listas según `what`.
  - `calculate_diversity`: aplica el cálculo a matrices, listas de `tximport` o objetos `SummarizedExperiment`, agregando por gen y devolviendo un `SummarizedExperiment` con assay `diversity` (S_q) o `hill` (D_q).

- Detección de diferencias y análisis estadístico:
  - `calculate_difference` y funciones en `difference_functions` calculan medias por grupo, diferencias (o log2 fold-changes), p-values y p-values ajustados. Soportan configuraciones típicas para análisis diferencial y están pensadas para integrarse con las funciones de diversidad.

- Wrappers y métodos auxiliares:
  - `calculate_method` gestiona la ejecución del método de diversidad por gen, formatea la salida y facilita evaluar múltiples `q` en una sola llamada.

- Inferencia y gestión de grupos de muestras:
  - `infer_sample_group` intenta deducir etiquetas de muestra (p. ej. Normal/Tumor, tipos TCGA) a partir de nombres de columnas o metadatos cuando no se proporcionan explícitamente.

- Visualizaciones integradas:
  - `plot_tsallis_q_curve`: muestra la mediana ± IQR de la entropía Tsallis a lo largo de diferentes `q` por grupo.
  - `plot_tsallis_violin_multq`: violines de Tsallis para múltiples `q` y grupos.
  - `plot_diversity_density`: densidades de diversidad por tipo de muestra.
  - `plot_mean_violin`: violines del valor medio por gen por tipo de muestra.
  - `plot_ma`: MA-plot para resultados diferenciales.
  - `plot_volcano`: volcano plot con etiquetado de los genes más significativos.
  - `plot_top_transcripts`: visualización/resumen de transcritos destacados por efecto.

- Utilidades internas y helpers:
  - `prepare_tsallis_long` y otros helpers transforman resultados para `ggplot2`, `tidyr` y `dplyr`, y facilitan la creación de gráficos y tablas resumen.

- Dataset de ejemplo incluido: `tcga_brca_luma_dataset`, incluido para reproducir ejemplos y pruebas.

Badges
------

- GitHub Actions (CI):

    [![R-CMD-check](https://github.com/gallardoalba/abyssEdge/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/gallardoalba/abyssEdge/actions)

Instalación
----------

- Instalar desde GitHub (recomendado durante desarrollo):

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("gallardoalba/abyssEdge")
```

Uso rápido
----------

```r
library(abyssEdge)

# Cargar el dataset incluido
data("tcga_brca_luma_dataset", package = "abyssEdge")

# Ejemplo: calcular diferencias entre grupos
res <- calculate_difference(tcga_brca_luma_dataset$counts,
                                                        group = tcga_brca_luma_dataset$group)
head(res)
```

Vignette y documentación
------------------------

- La viñeta está incluida en el paquete; puede verse con:

```r
browseVignettes("abyssEdge")
vignette("abyssEdge")
```

- La documentación completa puede consultarse en [inst/doc/abyssEdge.html](inst/doc/abyssEdge.html).

Pruebas y desarrollo
--------------------

- Ejecutar la batería de tests con `devtools::test()` o:

```sh
R -e 'devtools::test()'
```

- Para comprobaciones completas de paquete:

```sh
R CMD check --as-cran .
```

Citación y licencia
--------------------

- Ver el archivo de citación con `citation("abyssEdge")` o el fichero [inst/CITATION](inst/CITATION).
- Licencia: consulte [LICENSE](LICENSE).

Contribuir
----------

- Abra un issue o envíe un pull request en el repositorio GitHub.
- Mantenga el estilo de R/roxygen y ejecute `devtools::check()` antes de proponer cambios.

Soporte
------

Para preguntas o reportes de errores, abra un issue en el repositorio GitHub.

— El equipo de abyssEdge
````
