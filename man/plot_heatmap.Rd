% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions_results.R
\name{plot_heatmap}
\alias{plot_heatmap}
\title{Plot a heatmap}
\usage{
plot_heatmap(
  dep,
  type = c("contrast", "centered"),
  kmeans = FALSE,
  k = 6,
  col_limit = 6,
  indicate = NULL,
  clustering_distance = c("euclidean", "maximum", "manhattan", "canberra", "binary",
    "minkowski", "pearson", "spearman", "kendall", "gower"),
  row_font_size = 6,
  col_font_size = 10,
  plot = TRUE,
  ...
)
}
\arguments{
\item{dep}{SummarizedExperiment,
Data object for which differentially enriched proteins are annotated
(output from \code{\link{test_diff}()} and \code{\link{add_rejections}()}).}

\item{type}{'contrast' or 'centered',
The type of data scaling used for plotting.
Either the fold change ('contrast') or
the centered log2-intensity ('centered').}

\item{kmeans}{Logical(1),
Whether or not to perform k-means clustering.}

\item{k}{Integer(1),
Sets the number of k-means clusters.}

\item{col_limit}{Integer(1),
Sets the outer limits of the color scale.}

\item{indicate}{Character,
Sets additional annotation on the top of the heatmap
based on columns from the experimental design (colData).
Only applicable to type = 'centered'.}

\item{clustering_distance}{"euclidean", "maximum", "manhattan", "canberra",
"binary", "minkowski", "pearson", "spearman", "kendall" or "gower",
Function used to calculate clustering distance (for proteins and samples).
Based on \code{\link[ComplexHeatmap]{Heatmap}}
and \code{\link[cluster]{daisy}}.}

\item{row_font_size}{Integer(1),
Sets the size of row labels.}

\item{col_font_size}{Integer(1),
Sets the size of column labels.}

\item{plot}{Logical(1),
If \code{TRUE} (default) the heatmap is produced.
Otherwise (if \code{FALSE}), the data which the
heatmap is based on are returned.}

\item{...}{Additional arguments for Heatmap function as depicted in
\code{\link[ComplexHeatmap]{Heatmap}}}
}
\value{
A heatmap (generated by \code{\link[ComplexHeatmap]{Heatmap}})
}
\description{
\code{plot_heatmap} generates a heatmap of all significant proteins.
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

# Filter, normalize and impute missing values
filt <- filter_missval(se, thr = 0)
norm <- normalize_vsn(filt)
imputed <- impute(norm, fun = "MinProb", q = 0.01)

# Test for differentially expressed proteins
diff <- test_diff(imputed, "control", "Ctrl")
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

# Plot heatmap
plot_heatmap(dep)
plot_heatmap(dep, 'centered', kmeans = TRUE, k = 6, row_font_size = 3)
plot_heatmap(dep, 'contrast', col_limit = 10, row_font_size = 3)
}
