---
editor_options: 
  markdown: 
    wrap: 72
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# DEP2

<!-- badges: start -->

<!-- badges: end -->

DEP2 provides an comprehensive analysis workflow for mass spectrometry
based proteomics data, developed from the previous package DEP. This
package provided differential expression/enrichment analysis pipelines
for various data including protein-level quantity(e.g. proteingroup),
peptide-level quantity data and modification-specific proteomics
(quantities of modified peptides). Inherited from DEP, DEP2 provides
functions for data processing, hypothesis testing (via limma) and
visualization. The peptide to protein aggregation strategy is integrated
into workflows. To reduced the barrier in omics analysis, downstream
functional exploration are packaged as suites including functional
enrichment, timecourse cluster and protein-protein interaction network.
DEP2 also contains an easy-used shiny application designed under
modularization.

## Installation

You can install `DEP2` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("mildpiggy/DEP2")
```

## Run shiny application

Run the build-in shiny application though `run_app` function

``` r
library(DEP2)
DEP2::run_app()
```

## Vignettes

View the Vignettes for more information ...
