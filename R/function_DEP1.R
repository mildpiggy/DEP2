## here is the functions from DEP

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
#' data <- UbiLength
#' columns <- grep("LFQ.", colnames(data))
#'
#' # Get prefix
#' names <- colnames(data[, columns])
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
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' # Filter, normalize and impute missing values
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.01)
#'
#' # Test for differentially expressed proteins
#' diff <- test_diff(imputed, "control", "Ctrl")
#' dep <- add_rejections(diff, alpha = 0.05, lfc = 1)
#'
#' # Get a wide data.frame
#' wide <- get_df_wide(dep)
#' colnames(wide)
#' @export
get_df_wide <- function(se) {
  # Show error if inputs are not the required classes
  assert_that(inherits(se, "SummarizedExperiment"))

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
