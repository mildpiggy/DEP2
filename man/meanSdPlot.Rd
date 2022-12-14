% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions_QC.R
\name{meanSdPlot}
\alias{meanSdPlot}
\title{Plot row standard deviations versus row means}
\usage{
meanSdPlot(
  x,
  ranks = TRUE,
  xlab = ifelse(ranks, "rank(mean)", "mean"),
  ylab = "sd",
  pch,
  plot = TRUE,
  bins = 50,
  ...
)
}
\arguments{
\item{x}{SummarizedExperiment,
Data object.}

\item{ranks}{Logical,
Whether or not to plot the row means on the rank scale.}

\item{xlab}{Character,
x-axis label.}

\item{ylab}{Character,
y-axis label.}

\item{pch}{Ignored -
exists for backward compatibility.}

\item{plot}{Logical,
Whether or not to produce the plot.}

\item{bins}{Numeric vector,
Data object before normalization.}

\item{...}{Other arguments,
Passed to \code{\link[ggplot2:geom_hex]{stat_binhex}}.}
}
\value{
A scatter plot of row standard deviations
versus row means(generated by \code{\link[ggplot2:geom_hex]{stat_binhex}})
}
\description{
\code{meanSdPlot} generates a hexagonal heatmap
of the row standard deviations versus row means
from SummarizedExperiment objects.
See \code{\link[vsn]{meanSdPlot}}.
}
\examples{
# Load example
data <- UbiLength
data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Make SummarizedExperiment
columns <- grep("LFQ.", colnames(data_unique))
exp_design <- UbiLength_ExpDesign
se <- make_se(data_unique, columns, exp_design)

# Filter and normalize
filt <- filter_missval(se, thr = 0)
norm <- normalize_vsn(filt)

# Plot meanSdPlot
meanSdPlot(norm)
}
