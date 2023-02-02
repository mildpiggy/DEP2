% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEP_aggregation.R
\name{pe2se}
\alias{pe2se}
\title{Extract the proteins SummarizedExperiment object from a QFeatures container}
\usage{
pe2se(
  pe_aggregated,
  names = "Gene.names",
  ids = "smallestProteingroups",
  delim = ";"
)
}
\arguments{
\item{pe_aggregated}{A QFeatures object output from \code{aggregate_pe()}, which contain a "protein" quantative assay.}

\item{names}{The column of gene names, which serve as feature identifier and is transmitted to \link{make_unique}}

\item{ids}{The column of protein ID, transmitted to \link{make_unique}. aggregate_pe automatically generate a "smallestProteingroups" column
to store ids in proteingroups.}

\item{delim}{Character. The separator in names of ids.}
}
\value{
A SummarizedExperiment object.
}
\description{
This function accept a result from \link{aggregate_pe}(),
tidy features identifiers and return a SE object for following analyze
}
\examples{
\dontrun{
data <- Silicosis_peptide
ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = T)
pe_peptides <- make_pe_parse(Silicosis_peptide, columns = ecols, remove_prefix = T, log2transform = T,mode = "delim")
filt_pe <- filter_pe(pe_peptides, thr = 1,fraction = 0.4, filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" )
imp_pe <- QFeatures::addAssay(filt_pe, DEP2::impute(filt_pe[["peptideRaw"]], fun = "MinDet"), name = "peptideImp")
norm_pe <- DEP2:::normalize_pe(imp_pe,method = "quantiles", i = "peptideImp", name = "peptideNorm")
protein_pe <- DEP2::aggregate_pe(norm_pe, fcol = "Proteins", peptide_assay_name = "peptideNorm")
class(protein_pe)
se <- pe2se(protein_pe)
class(se)
}
}