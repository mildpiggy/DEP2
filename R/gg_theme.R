#' DEP ggplot theme 1
#'
#' \code{theme_DEP1} is the default ggplot theme used for plotting
#' in \code{\link{DEP2}} with horizontal x-axis labels.
#' @return ggplot theme
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != '+' & data$Potential.contaminant != '+',]
#' data_unique <- make_unique(data, 'Gene.names', 'Protein.IDs', delim = ';')
#'
#' columns <- grep('LFQ.', colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' plot_frequency(filt) # uses theme_DEP1() style
#'
#' @export
theme_DEP1 <- function() {
    # Use theme_bw() as default
    basesize <- 12
    theme <- ggplot2::theme_bw(base_size = basesize)

    # Change plot title appearance
    theme$plot.title$face <- "bold"
    theme$plot.title$size <- basesize + 2
    theme$plot.title$hjust <- 0.5

    # Change axis title appearance
    theme$axis.title.x$size <- basesize + 2

    theme$axis.title.y$size <- basesize + 2

    # Change axis text appearance
    theme$axis.text$size <- basesize
    theme$axis.text$colour <- "black"

    # Change legend title appearance
    theme$legend.title$size <- basesize + 2

    # Change legend text appearance
    theme$legend.text$size <- basesize

    # Change strip text (facet headers) appearance
    theme$strip.text$face <- "bold"
    theme$strip.text$size <- basesize + 2
    theme$strip.text$colour <- "black"

    return(theme)
}

#' DEP ggplot theme 2
#'
#' \code{theme_DEP2} is the ggplot theme used for plotting
#' in \code{\link{DEP2}} with vertical x-axis labels.
#' @return ggplot theme
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != '+' & data$Potential.contaminant != '+',]
#' data_unique <- make_unique(data, 'Gene.names', 'Protein.IDs', delim = ';')
#'
#' columns <- grep('LFQ.', colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' plot_numbers(filt) # uses theme_DEP2() style
#'
#' @export
#' @export
theme_DEP2 <- function() {
    # Get vertical x-axis labels
    theme <- theme_DEP1()
    theme$axis.text.x$angle <- 90
    theme$axis.text.x$hjust <- 1
    theme$axis.text.x$vjust <- 0.5
    return(theme)
}
