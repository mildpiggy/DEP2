% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_DEP1.R
\name{plot_frequency}
\alias{plot_frequency}
\title{Plot protein overlap between samples}
\usage{
plot_frequency(se, plot = TRUE)
}
\arguments{
\item{se}{SummarizedExperiment,
Data object for which to plot observation frequency.}

\item{plot}{Logical(1),
If \code{TRUE} (default) the barplot is produced.
Otherwise (if \code{FALSE}), the data which the
barplot is based on are returned.}
}
\value{
Barplot of overlap of protein identifications
between samples (generated by \code{\link[ggplot2]{ggplot}})
}
\description{
\code{plot_frequency} generates a barplot
of the protein overlap between samples
}
\examples{
# Load example
data <- Silicosis_pg
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
# Make SummarizedExperiment
ecols <- grep("LFQ.", colnames(data_unique))
se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")

# Filter and plot frequency
filt <- filter_se(se, thr = 0, fraction = 0.3)

plot_frequency(filt)
}