# TSENAT

R-Paket zum Erkennen und Analysieren von
Expressions-/Transkriptionsunterschieden und zur Berechnung von
Diversitätsmaßen (enthält Tsallis-Funktionen, statistische Tests und
Visualisierungswerkzeuge).

## Herkunft und Zuordnung

TSENAT baut auf wesentlichen Teilen des Codes des Projekts
SplicingFactory auf (Anerkennung der ursprünglichen Autoren). Der Code
wurde angepasst und um zusätzliche Werkzeuge, Fehlerbehebungen und
Visualisierungsfunktionen für die Analyse der Transkriptvielfalt
erweitert.

## Tsallis-Theorie (kurz)

Die Tsallis-Entropie ist eine Verallgemeinerung der Shannon-Entropie und
wird definiert als S_q = (1 - sum p^q) / (q - 1) für eine
Wahrscheinlichkeitsverteilung p. Im Grenzfall q -\> 1 ergibt sich die
Shannon-Entropie. In TSENAT wird die Tsallis-Entropie sowie die
zugehörigen Hill-Zahlen (D_q) pro Gen berechnet, um Isoform-Vielfalt zu
quantifizieren. Der Parameter `q` steuert die Empfindlichkeit gegenüber
seltenen vs. häufigen Isoformen (q \< 1 betont Seltene, q \> 1 betont
Häufige).

## Integrierte Funktionalitäten

- Berechnung von Tsallis-Entropie und Diversität:
  - `calculate_tsallis_entropy`: berechnet S_q und/oder D_q für einen
    numerischen Vektor (unterstützt Normalisierung, mehrere `q`-Werte
    und den Grenzfall q→1). Rückgabe als numerischer Vektor oder Liste
    je nach `what`.
  - `calculate_diversity`: führt die Berechnung über Transkripte/Genera
    für Matrizen, `tximport`-Listen oder `SummarizedExperiment`-Objekte
    aus und gibt ein `SummarizedExperiment` mit dem Assay `diversity`
    (S_q) oder `hill` (D_q) zurück.
- Differenz- und Statistikfunktionen:
  - `calculate_difference` und Hilfsfunktionen in `difference_functions`
    berechnen Gruppenmittelwerte, Unterschiede (oder log2 Fold-Changes),
    p-Werte und adjustierte p-Werte. Diese Funktionen sind so gestaltet,
    dass sie mit Diversitätszusammenfassungen sowie mit
    Expressionsmatrizen zusammenarbeiten.
- Wrapper und Hilfsfunktionen:
  - `calculate_method` bietet einen Wrapper, um die gewählte
    Diversitätsmethode pro Gen auszuführen, Ausgaben zu formatieren und
    mehrere `q`-Werte in einem Durchlauf auszuwerten.
- Gruppenerkennung und Metadaten-Hilfen:
  - `infer_sample_group` versucht, Stichprobenklassen (z. B.
    Normal/Tumor oder TCGA-Codes) aus Spaltennamen oder Metadaten
    abzuleiten, wenn keine expliziten Gruppenlabels vorhanden sind.
- Visualisierung und Plots:
  - `plot_tsallis_q_curve`: Median ± IQR der Tsallis-Entropie über
    q-Werte getrennt nach Gruppen.
  - `plot_tsallis_violin_multq`: Violin-Plots der Tsallis-Entropie für
    mehrere q-Werte und Gruppen.
  - `plot_diversity_density`: Dichteplots der Diversität nach
    Stichprobentyp.
  - `plot_mean_violin`: Violin-Plot des pro-Gen-Mittelwerts der
    Diversität nach Stichprobentyp.
  - `plot_ma`: MA-Plot für differentielle Ergebnisse.
  - `plot_volcano`: Volcano-Plot mit Beschriftung der wichtigsten Gene.
  - `plot_top_transcripts`: Zusammenfassung/Visualisierung der
    wichtigsten Transkripte.
- Interne Helfer:
  - Hilfsfunktionen wie `prepare_tsallis_long` formatieren Ergebnisse
    für `ggplot2` und arbeiten mit `tidyr`/`dplyr`, um sofort
    darstellbare Data-Frames zu erzeugen.
- Beispieldatensatz: `tcga_brca_luma_dataset`, enthalten für
  Vignetten-Beispiele und Tests.

## Installation

Installation von GitHub (während der Entwicklung empfohlen):

``` r

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/TSENAT")
```

## Schnellstart

``` r

library(TSENAT)
data("tcga_brca_luma_dataset", package = "TSENAT")

res <- calculate_difference(tcga_brca_luma_dataset$counts,
                            group = tcga_brca_luma_dataset$group)
head(res)
```

## Vignette und Dokumentation

Die Vignette ist im Paket enthalten; anzeigen mit:

``` r

browseVignettes("TSENAT")
vignette("TSENAT")
```

## Tests und Entwicklung

Tests ausführen mit
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
oder:

``` sh
R -e 'devtools::test()'
```

Paket prüfen:

``` sh
R CMD check --as-cran .
```

## Zitation und Lizenz

Siehe `citation("TSENAT")` oder
[inst/CITATION](https://gallardoalba.github.io/TSENAT/inst/CITATION).
Lizenz in [LICENSE](https://gallardoalba.github.io/TSENAT/LICENSE).

## Beitragen

Öffnen Sie ein Issue oder senden Sie einen Pull Request auf GitHub.
Bitte führen Sie
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
vor dem Einreichen aus.

## Veröffentlichen auf GitHub (kurze Schritte)

1.  Remote-Repo erstellen und pushen:

``` sh
git init
git add .
git commit -m "Initial commit: add package source and README"
gh repo create YOUR_GITHUB_USERNAME/TSENAT --public --source=. --remote=origin
git push -u origin main
```

2.  CI hinzufügen: die enthaltene GitHub Actions Workflow-Datei nutzen
    oder `usethis::use_github_action_check_standard()` aufrufen.

## Support

Fehler bitte im GitHub-Repository melden.

— Das TSENAT-Team
