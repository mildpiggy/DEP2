## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align="center")

## ----run,eval=FALSE-----------------------------------------------------------
#  DEP2::run_app()

## ----warning=FALSE,message=FALSE----------------------------------------------
## Check required packages for enrichment analysis
DEP2::check_enrichment_depends()

## Check required packages for PPI
DEP2::check_PPI_depends()

## ----eval=FALSE---------------------------------------------------------------
#  ## Install all depend packages at once by setting 'dependencies = TRUE'
#  devtools::install_github("mildpiggy/DEP2", dependencies = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  ## Check organism annotation package of a species
#  ## Parameter organism could one of "Chicken", "Human", "Mouse", "Rat", etc.
#  DEP2::check_organismDB_depends(organism = "Human")

## ----eval=FALSE---------------------------------------------------------------
#  DEP2::check_PPI_depends()

