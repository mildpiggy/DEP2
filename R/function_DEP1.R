## Here is the functions entirely inherited from DEP

#' Obtain the longest common prefix
#'
#' \code{get_prefix} returns the longest common prefix
#' of the supplied words.
#'
#' @param words Character vector,
#' A list of words.
#' @return A character vector containing the prefix.
#' @examples
#' # Load example
#' data <- Silicosis_pg
#' names <- grep("LFQ.", colnames(data), value = TRUE)
#'
#' # Get prefix
#' get_prefix(names)
#' @export
get_prefix <- function(words) {
  # Show error if input is not the required class
  assertthat::assert_that(is.character(words))

  # Show error if 'words' contains 1 or less elements
  if(length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  # Show error if 'words' contains NA
  if(any(is.na(words))) {
    stop("'words' contains NAs")
  }

  # Truncate words to smallest name
  minlen <- min(nchar(words))
  truncated <- substr(words, 1, minlen)

  # Show error if one of the elements is shorter than one character
  if(minlen < 1) {
    stop("At least one of the elements is too short")
  }

  # Get identifical characters
  mat <- data.frame(strsplit(truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 1)

  # Obtain the longest common prefix
  prefix <- as.logical(cumprod(identical))
  paste(mat[prefix, 1], collapse = "")
}

#' Obtain the longest common suffix
#'
#' \code{get_suffix} returns the longest common suffix
#' of the supplied words.
#'
#' @param words Character vector,
#' A list of words.
#' @return A character vector containing the suffix
#' @examples
#' # Get suffix
#' names <- c("xyz_rep", "abc_rep")
#' get_suffix(names)
#' @export
get_suffix <- function(words) {
  # Show error if input is not the required class
  assertthat::assert_that(is.character(words))

  # Show error if 'words' contains 1 or less elements
  if(length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  # Show error if 'words' contains NA
  if(any(is.na(words))) {
    stop("'words' contains NAs")
  }

  # Truncate words to smallest name
  minlen <- min(nchar(words))
  truncated <- substr(words, nchar(words) - minlen + 1, nchar(words))

  # Show error if one of the elements is shorter than one character
  if(minlen < 1) {
    stop("At least one of the elements is too short")
  }

  # Reverse characters wihtin word
  rev_string <- function(str) {
    paste(rev(strsplit(str, "")[[1]]), collapse = "")
  }
  rev_truncated <- vapply(truncated, rev_string, character(1))

  # Get identifical characters
  mat <- data.frame(strsplit(rev_truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 1)

  # Obtain the longest common prefix
  prefix <- as.logical(cumprod(identical))
  rev_string(paste(mat[prefix, 1], collapse = ""))
}

# Short internal function to delete the longest common prefix
delete_prefix <- function(words) {
  # Get prefix
  prefix <- get_prefix(words)
  # Delete prefix from words
  gsub(paste0("^", prefix), "", words)
}

# Short internal function to delete the longest common suffix
delete_suffix <- function(words) {
  # Get prefix
  suffix <- get_suffix(words)
  # Delete prefix from words
  gsub(paste0(suffix, "$"), "", words)
}



#' Generate a wide data.frame from a SummarizedExperiment
#'
#' \code{get_df_wide} generate a wide data.frame from a SummarizedExperiment.
#'
#' @param se SummarizedExperiment or DEGdata object.
#' @return A data.frame object
#' containing all data in a wide format(combine the assay and rowData),
#' where each row represents a protein/gene.
#' @examples
#' # Load example
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")
#' # Filter and normalize
#' filt <- filter_se(se, thr = 0, fraction = 0.4, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinDet")
#' diff <- DEP2::test_diff(imputed,type = "control", control = "PBS")
#' dep <- add_rejections(diff)
#' get_df_wide(dep)
#'
#' @export
get_df_wide <- function(se) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # Show error if inputs do not contain required columns
  if (!"name" %in% colnames(rowData(se, use.names = FALSE))) {
    stop("'name' column is not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }

  # Extract row data
  row_data <- rowData(se, use.names = FALSE) %>%
    data.frame()
  # Extract assay data
  assay_data <- assay(se) %>%
    data.frame() %>%
    rownames_to_column()
  colnames(assay_data)[1] <- "name"

  # Merge row and assay data into a wide data.frame
  wide <- full_join(assay_data, row_data, by = "name")

  return(wide)
}


#' Plot protein overlap between samples
#'
#' \code{plot_frequency} generates a barplot
#' of the protein overlap between samples
#'
#' @param se SummarizedExperiment,
#' Data object for which to plot observation frequency.
#' @param plot Logical(1),
#' If \code{TRUE} (default) the barplot is produced.
#' Otherwise (if \code{FALSE}), the data which the
#' barplot is based on are returned.
#' @return Barplot of overlap of protein identifications
#' between samples (generated by \code{\link[ggplot2]{ggplot}})
#' @examples
#' # Load example
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")
#'
#' # Filter and plot frequency
#' filt <- filter_se(se, thr = 0, fraction = 0.3)
#'
#' plot_frequency(filt)
#' @export
plot_frequency <- function(se, plot = TRUE) {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot),
                          length(plot) == 1)

  # Make a binary long data.frame (1 = valid value, 0 = missing value)
  df <- assay(se) %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(ID, bin, -rowname) %>%
    mutate(bin = ifelse(is.na(bin), 0, 1))
  # Identify the number of experiments a protein was observed
  stat <- df %>%
    group_by(rowname) %>%
    summarize(sum = sum(bin))
  # Get the frequency of the number of experiments proteins
  # were observerd and plot these numbers
  table <- table(stat$sum) %>% data.frame()
  p <- ggplot(table, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() +
    scale_fill_grey(start = 0.8, end = 0.2) +
    labs(title = "Protein identifications overlap",
         x = "Identified in number of samples",
         y = "Number of proteins") +
    theme_DEP2() +
    theme(legend.position="none")
  if(plot) {
    return(p)
  } else {
    df <- as.data.frame(table)
    colnames(df) <- c("samples", "proteins")
    return(df)
  }
}


#' Plot frequency of significant conditions per protein
#'
#' \code{plot_cond_freq} generates a histogram of the number of significant conditions per protein.
#'
#' @param dep SummarizedExperiment,
#' Data object for which differentially enriched proteins are annotated
#' (output from \code{\link{test_diff}()} and \code{\link{add_rejections}()}).
#' @param plot Logical(1),
#' If \code{TRUE} (default) the histogram is produced.
#' Otherwise (if \code{FALSE}), the data which the
#' histogram is based on are returned.
#' @return A histogram (generated by \code{\link[ggplot2]{ggplot}})
#' @examples
#' # Load example
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")
#'
#' # Analysis
#' filt <- filter_se(se, thr = 0, fraction = 0.3, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinDet")
#' diff <- test_diff(imputed,type = "control", control = "PBS")
#' dep <- add_rejections(diff)
#'
#' # Plot frequency of significant conditions
#' plot_cond_freq(dep)
#' @export
plot_cond_freq <- function(dep, plot = TRUE) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(plot),
                          length(plot) == 1)

  # Check for significance columns
  row_data <- rowData(dep, use.names = FALSE)
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }

  # Check for significant column
  if(!"significant" %in% colnames(row_data)) {
    stop("'significant' column is not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required column",
         call. = FALSE)
  }

  # Filter for significant proteins
  significant <- dep[row_data$significant, ]

  # Get significant columns and count significant conditions per protein
  row_data <- rowData(significant, use.names = FALSE) %>% data.frame()
  cols <- grep("_significant", colnames(row_data))
  df <- row_data %>%
    select(name, ID, cols) %>%
    gather(condition, significant, -c(name, ID)) %>%
    mutate(val = ifelse(significant, 1, 0))
  stat <- df %>%
    group_by(name) %>%
    summarize(sum = sum(val), n = n())
  # Get a counts table
  table <- table(stat$sum) %>% data.frame()

  # Plot the count table as a bar graph
  p <- ggplot(table, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() +
    labs(x = "Number of significant conditions",
         y = "Number of proteins",
         title = "Frequency of significant conditions per protein") +
    DEP::theme_DEP1() +
    theme(legend.position = "none")
  if(plot) {
    return(p)
  } else {
    df <- as.data.frame(table)
    colnames(df) <- c("conditions", "proteins")
    return(df)
  }
}

#' Plot conditions overlap
#'
#' \code{plot_cond_overlap} generates a histogram of
#' the number of proteins per condition or overlapping conditions.
#'
#' @param dep SummarizedExperiment,
#' Data object for which differentially enriched proteins are annotated
#' (output from \code{\link{test_diff}()} and \code{\link{add_rejections}()}).
#' @param plot Logical(1),
#' If \code{TRUE} (default) the barplot is produced.
#' Otherwise (if \code{FALSE}), the data which the
#' barplot is based on are returned.
#' @return A histogram (generated by \code{\link[ggplot2]{ggplot}})
#' @examples
#' # Load example
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")
#'
#' # Analysis
#' filt <- filter_se(se, thr = 0, fraction = 0.3, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinDet")
#' diff <- test_diff(imputed,type = "control", control = "PBS")
#' dep <- add_rejections(diff)
#'
#' # Plot condition overlap
#' plot_cond_overlap(dep)
#' @export
plot_cond_overlap <- function(dep, plot = TRUE) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(plot),
                          length(plot) == 1)

  # Check for significance columns
  row_data <- rowData(dep, use.names = FALSE)
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }

  # Check for significant column
  if(!"significant" %in% colnames(row_data)) {
    stop("'significant' column is not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required column",
         call. = FALSE)
  }

  # Filter for significant proteins
  significant <- dep[row_data$significant, ]

  # Get significant columns
  row_data <- rowData(significant, use.names = FALSE) %>% data.frame()
  cols <- grep("_significant", colnames(row_data))
  colnames(row_data)[cols] <- gsub("_significant", "", colnames(row_data)[cols])

  # Rename column names
  row_data_renamed <- row_data
  colnames(row_data_renamed)[cols] <- LETTERS[seq(to = length(cols))]
  legend <- data.frame(symbol = colnames(row_data_renamed)[cols],
                       contrast = colnames(row_data)[cols])

  # Get co-occurence matrix
  df <- select(row_data_renamed, name, ID, cols)
  counts <- table(df[,3:(length(cols)+2)]) %>%
    as.data.frame() %>%
    filter(Freq > 0)

  # Parse condition names
  mat <- counts[,seq_len(ncol(counts)-1)] == TRUE
  counts$conditions <-
    apply(mat, 1, function(x) {
      paste0(colnames(counts[,seq_len(ncol(counts)-1)])[x], collapse = " ")
    })

  # Sort on number of conditions
  counts$n_con <- apply(mat, 1, function(x) length(which(x)))
  counts <- counts %>%
    arrange(n_con, conditions)
  counts$conditions <- readr::parse_factor(counts$conditions, levels = counts$conditions)

  if(nrow(counts) <= 40 ) {
    labelsize = 12
  } else {
    labelsize = 12 / (nrow(counts) / 40)
  }

  # Plot conditions overlap histogram
  p1 <- ggplot(counts, aes(x = conditions, y = Freq)) +
    geom_col(fill = "black") +
    labs(title = "Overlap between conditions",
         x = "Conditions",
         y = "Number of Proteins") +
    theme_DEP2() +
    theme(axis.text=element_text(size = labelsize))

  # Legend table
  ttheme <- gridExtra::ttheme_minimal(
    core=list(fg_params=list(hjust=0, x=0.1)))
  p2 <- gridExtra::tableGrob(legend, theme = ttheme, rows = NULL, cols = NULL)
  if(!plot) {
    df <- counts %>%
      select(conditions, Freq, n_con)
    colnames(df) <- c("conditions", "proteins", "number_of_conditions")
    list <- list(counts = df, legend = legend)
    return(list)
  } else {
    gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(0.8, 0.2))
  }
}

#' Plot frequency of significant conditions per protein
#' and the overlap in proteins between conditions
#'
#' \code{plot_cond} generates a histogram of
#' the number of proteins per condition and stacks for overlapping conditions.
#' The function is inherited from DEP package.
#'
#' @param dep SummarizedExperiment,
#' Data object for which differentially enriched proteins are annotated
#' (output from \code{\link{test_diff}()} and \code{\link{add_rejections}()}).
#' @param plot Logical(1),
#' If \code{TRUE} (default) the barplot is produced.
#' Otherwise (if \code{FALSE}), the data which the
#' barplot is based on are returned.
#' @return A histogram (generated by \code{\link[ggplot2]{ggplot}})
#' @examples
#' # Load example
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")
#'
#' # Analysis
#' filt <- filter_se(se, thr = 0, fraction = 0.3, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinDet")
#' diff <- test_diff(imputed,type = "control", control = "PBS")
#' dep <- add_rejections(diff)
#'
#' # Plot histogram with overlaps
#' plot_cond(dep)
#' @export
plot_cond <- function(dep, plot = TRUE) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(plot),
                          length(plot) == 1)

  # Check for significance columns
  row_data <- rowData(dep, use.names = FALSE)
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }

  # Check for significant column
  if(!"significant" %in% colnames(row_data)) {
    stop("'significant' column is not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required column",
         call. = FALSE)
  }

  # Filter for significant proteins
  significant <- dep[row_data$significant, ]

  # Get significant columns
  row_data <- rowData(significant, use.names = FALSE) %>% data.frame()
  cols <- grep("_significant", colnames(row_data))
  colnames(row_data)[cols] <- gsub("_significant", "", colnames(row_data)[cols])

  # Rename column names
  row_data_renamed <- row_data
  colnames(row_data_renamed)[cols] <- LETTERS[seq(to = length(cols))]
  legend <- data.frame(symbol = colnames(row_data_renamed)[cols], names = colnames(row_data)[cols])

  # Get co-occurence matrix
  df <- select(row_data_renamed, name, ID, cols)
  counts <- table(df[,3:(length(cols)+2)]) %>%
    as.data.frame() %>%
    filter(Freq > 0)

  # Parse condition names
  mat <- counts[,seq_len(ncol(counts)-1)] == TRUE
  counts$conditions <-
    apply(mat, 1, function(x) {
      paste0(colnames(counts[,seq_len(ncol(counts)-1)])[x], collapse = " ")
    })

  # Sort on number of conditions
  counts$n_con <- apply(mat, 1, function(x) length(which(x)))
  counts <- counts %>%
    arrange(n_con, conditions)
  counts$conditions <- readr::parse_factor(counts$conditions, levels = counts$conditions)

  counts <- mutate(counts, ID = paste0(conditions, ": ", Freq))

  if(nrow(counts) <= 40 ) {
    labelsize = 12
  } else {
    labelsize = 12 / (nrow(counts) / 40)
  }

  # Plot conditions overlap histogram
  p1 <- ggplot(counts, aes(x = n_con, y = Freq, fill = conditions)) +
    geom_col() +
    geom_text(aes(label = ID), size = 3, hjust = 0.5,
              position = position_stack(vjust = 0.5)) +
    labs(x = "Number of significant conditions",
         y = "Number of proteins",
         title = "Overlap between conditions") +
    DEP::theme_DEP1()

  # Legend table
  ttheme <- gridExtra::ttheme_minimal(
    core=list(fg_params=list(hjust=0, x=0.1)))
  p2 <- gridExtra::tableGrob(legend, theme = ttheme, rows = NULL, cols = NULL)
  if(!plot) {
    df <- counts %>%
      select(conditions, Freq, n_con)
    colnames(df) <- c("conditions", "proteins", "number_of_conditions")
    list <- list(counts = df, legend = legend)
    return(list)
  } else{
    gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(0.8, 0.2))
  }
}


#' DEP ggplot theme 1
#'
#' \code{theme_DEP1} is the default ggplot theme inherited  used for plotting
#' in DEP package with horizontal x-axis labels.
#' @return ggplot theme
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
#' \code{theme_DEP2} is the ggplot theme from DEP package used for plotting
#' in DEP package with vertical x-axis labels.
#' @return ggplot theme
#' @export
theme_DEP2 <- function() {
  # Get vertical x-axis labels
  theme <- theme_DEP1()
  theme$axis.text.x$angle <- 90
  theme$axis.text.x$hjust <- 1
  theme$axis.text.x$vjust <- 0.5
  return(theme)
}
