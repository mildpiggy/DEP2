% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gg_theme.R
\name{theme_DEP1}
\alias{theme_DEP1}
\title{DEP ggplot theme 1}
\usage{
theme_DEP1()
}
\value{
ggplot theme
}
\description{
\code{theme_DEP1} is the default ggplot theme used for plotting
in \code{\link{DEP2}} with horizontal x-axis labels.
}
\examples{
data <- UbiLength
data <- data[data$Reverse != '+' & data$Potential.contaminant != '+',]
data_unique <- make_unique(data, 'Gene.names', 'Protein.IDs', delim = ';')

columns <- grep('LFQ.', colnames(data_unique))
exp_design <- UbiLength_ExpDesign
se <- make_se(data_unique, columns, exp_design)

filt <- filter_missval(se, thr = 0)
plot_frequency(filt) # uses theme_DEP1() style

}
