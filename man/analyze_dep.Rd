% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_functions.R
\name{analyze_dep}
\alias{analyze_dep}
\title{Differential expression analysis}
\usage{
analyze_dep(
  se,
  type = c("all", "control", "manual"),
  control = NULL,
  alpha = 0.05,
  lfc = 1,
  test = NULL,
  design_formula = formula(~0 + condition)
)
}
\arguments{
\item{se}{SummarizedExperiment,
Proteomics data with unique names and identifiers
annotated in 'name' and 'ID' columns.
Additionally, the colData should contain sample annotation including
'label', 'condition' and 'replicate' columns.
The appropriate columns and objects can be generated
using \code{\link{make_se}} or \code{\link{make_se_parse}}.}

\item{type}{"all", "control" or "manual",
The type of contrasts that will be tested.
This can be all possible pairwise comparisons ("all"),
limited to the comparisons versus the control ("control"), or
manually defined contrasts ("manual").}

\item{control}{Character(1),
The condition to which contrasts are generated
(a control condition would be most appropriate).}

\item{alpha}{Numeric(1),
Sets the threshold for the adjusted P value.}

\item{lfc}{Numeric(1),
Sets the threshold for the log2 fold change.}

\item{test}{Character,
The contrasts that will be tested if type = "manual".
These should be formatted as "SampleA_vs_SampleB" or
c("SampleA_vs_SampleC", "SampleB_vs_SampleC").}

\item{design_formula}{Formula,
Used to create the design matrix.}
}
\value{
A SummarizedExperiment object
containing FDR estimates of differential expression and
logical columns indicating significant proteins.
}
\description{
\code{analyze_dep} tests for differential expression of proteins
based on protein-wise linear models and empirical Bayes
statistics using \pkg{limma}.
}
\examples{
# Load datasets
data <- UbiLength
exp_design <- UbiLength_ExpDesign

# Import and process data
se <- import_MaxQuant(data, exp_design)
processed <- process(se)

# Differential protein expression analysis
dep <- analyze_dep(processed, "control", "Ctrl")
dep <- analyze_dep(processed, "control", "Ctrl",
    alpha = 0.01, lfc = log2(1.5))
dep <- analyze_dep(processed, "manual", test = c("Ubi6_vs_Ubi4"))
}
