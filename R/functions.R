#' Make unique names
#'
#' \code{make_unique} generates unique identifiers
#' for a proteomics dataset based on "name" and "id" columns.
#'
#' @param proteins Data.frame,
#' Protein table for which unique names will be created.
#' @param names Character(1),
#' Name of the column containing feature names.
#' @param ids Character(1),
#' Name of the column containing feature IDs.
#' @param delim Character(1),
#' Sets the delimiter separating the feature names within one protein group.
#' @return A data.frame with the additional variables
#' "name" and "ID" containing unique names and identifiers, respectively.
#' @examples
#' # Load example
#' data <- UbiLength
#'
#' # Check colnames and pick the appropriate columns
#' colnames(data)
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#' @export
make_unique <- function(proteins, names, ids, delim = ";") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins),
    is.character(names),
    length(names) == 1,
    is.character(ids),
    length(ids) == 1,
    is.character(delim),
    length(delim) == 1)

  col_names <- colnames(proteins)
  # Show error if inputs do not contain required columns
  if(!names %in% col_names) {
    stop("'", names, "' is not a column in '",
         deparse(substitute(proteins)), "'",
         call. = FALSE)
  }
  if(!ids %in% col_names) {
    stop("'", ids, "' is not a column in '",
         deparse(substitute(proteins)), "'",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins))
    proteins <- as.data.frame(proteins)

  # Select the name and id columns, and check for NAs
  double_NAs <- apply(proteins[,c(names, ids)], 1, function(x) all(is.na(x)))
  if(any(double_NAs)) {
    stop("NAs in both the 'names' and 'ids' columns")
  }

  # Take the first identifier per row and make unique names.
  # If there is no name, the ID will be taken.
  proteins_unique <- proteins %>%
    mutate(name = gsub(paste0(delim, ".*"), "", get(names)),
      ID = gsub(paste0(delim, ".*"), "", get(ids)),
      name = make.unique(ifelse(name == "" | is.na(name), ID, name)))
  return(proteins_unique)
}

#' Data.frame to SummarizedExperiment object
#' conversion using an experimental design
#'
#' \code{make_se} creates a SummarizedExperiment object
#' based on two data.frames: the protein table and experimental design.
#'
#' @param proteins_unique Data.frame,
#' Protein table with unique names annotated in the 'name' column
#' (output from \code{\link{make_unique}()}).
#' @param columns Integer vector,
#' Column numbers indicating the columns containing the assay data.
#' @param expdesign Data.frame,
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' See \code{\link{UbiLength_ExpDesign}} for an example experimental design.
#' @return A SummarizedExperiment object
#' with log2-transformed values.
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
#' @export
make_se <- function(proteins_unique, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.data.frame(expdesign))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if(tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)

  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)

  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID

  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }

  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]

  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name

  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                                  colData = expdesign,
                                  rowData = row_data)
  return(se)
}

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

#' Data.frame to SummarizedExperiment object
#' conversion using parsing from column names
#'
#' \code{make_se_parse} creates a SummarizedExperiment object
#' based on a single data.frame.
#'
#' @param proteins_unique Data.frame,
#' Protein table with unique names annotated in the 'name' column
#' (output from \code{\link{make_unique}()}).
#' @param columns Integer vector,
#' Column numbers indicating the columns containing the assay data.
#' @param mode "char" or "delim",
#' The mode of parsing the column headers.
#' "char" will parse the last number of characters as replicate number
#' and requires the 'chars' parameter.
#' "delim" will parse on the separator and requires the 'sep' parameter.
#' @param chars Numeric(1),
#' The number of characters to take at the end of the column headers
#' as replicate number (only for mode == "char").
#' @param sep Character(1),
#' The separator used to parse the column header
#' (only for mode == "delim").
#' @return A SummarizedExperiment object
#' with log2-transformed values.
#' @examples
#' # Load example
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' columns <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, columns, mode = "char", chars = 1)
#' se <- make_se_parse(data_unique, columns, mode = "delim", sep = "_")
#' @export
make_se_parse <- function(proteins_unique, columns,
                          mode = c("char", "delim"), chars = 1, sep = "_") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.character(mode),
                          is.numeric(chars),
                          length(chars) == 1,
                          is.character(sep),
                          length(sep) == 1)

  # Show error if inputs do not contain required columns
  mode <- match.arg(mode)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)

  # Select the assay values
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  colnames(raw) <- delete_prefix(colnames(raw)) %>% make.names()

  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name

  # Generate the colData
  if(mode == "char") {
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      mutate(condition = substr(label, 1, nchar(label) - chars),
             replicate = substr(label, nchar(label) + 1 - chars, nchar(label))) %>%
      unite(ID, condition, replicate, remove = FALSE)
  }
  if(mode == "delim") {
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      separate(label, c("condition", "replicate"), sep = sep,
               remove = FALSE, extra = "merge") %>%
      unite(ID, condition, replicate, remove = FALSE)
  }
  rownames(col_data) <- col_data$ID
  colnames(raw)[match(col_data$label, colnames(raw))] <- col_data$ID
  raw <- raw[, !is.na(colnames(raw))]


  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                             colData = col_data,
                             rowData = row_data)
  return(se)
}


#' Filter on missing values
#'
#' \code{filter_missval} filters a proteomics dataset based on missing values.
#' The dataset is filtered for proteins that have a maximum of
#' 'thr' missing values in at least one condition.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}).
#' @param thr Integer(1),
#' Sets the threshold for the allowed number of missing values
#' in at least one condition.
#' @return A filtered SummarizedExperiment object.
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
#' # Filter
#' stringent_filter <- filter_missval(se, thr = 0)
#' less_stringent_filter <- filter_missval(se, thr = 1)
#' @export
filter_missval <- function(se, thr = 0) {
  # Show error if inputs are not the required classes
  if(is.integer(thr)) thr <- as.numeric(thr)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(thr),
                          length(thr) == 1)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  max_repl <- max(colData(se)$replicate)
  if(thr < 0 | thr > max_repl) {
    stop("invalid filter threshold applied",
         "\nRun filter_missval() with a threshold ranging from 0 to ",
         max_repl)
  }

  # Make assay values binary (1 = valid value)
  bin_data <- assay(se)
  idx <- is.na(assay(se))
  bin_data[!idx] <- 1
  bin_data[idx] <- 0

  # Filter se on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(ID, value, -rowname) %>%
    left_join(., data.frame(colData(se)), by = "ID") %>%
    group_by(rowname, condition) %>%
    summarize(miss_val = n() - sum(value)) %>%
    filter(miss_val <= thr) %>%
    spread(condition, miss_val)
  se_fltrd <- se[keep$rowname, ]
  return(se_fltrd)
}

#' Filter proteins based on missing values
#'
#' \code{filter_proteins} filters a proteomic dataset based on missing values.
#' Different types of filtering can be applied, which range from only keeping
#' proteins without missing values to keeping proteins with a certain percent
#' valid values in all samples or keeping proteins that are complete
#' in at least one condition.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}).
#' @param type "complete", "condition" or "fraction",
#' Sets the type of filtering applied. "complete" will only keep
#' proteins with valid values in all samples. "condition" will keep
#' proteins that have a maximum of 'thr' missing values in at least
#' one condition. "fraction" will keep proteins that have a certain
#' fraction of valid values in all samples.
#' @param thr Integer(1),
#' Sets the threshold for the allowed number of missing values
#' in at least one condition if type = "condition".
#' @param min Numeric(1),
#' Sets the threshold for the minimum fraction of valid values
#' allowed for any protein if type = "fraction".
#' @return A filtered SummarizedExperiment object.
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
#' # Filter
#' stringent_filter <- filter_proteins(se, type = "complete")
#' less_stringent_filter <- filter_proteins(se, type = "condition", thr = 0)
#' @export
filter_proteins <- function(se, type = c("complete", "condition", "fraction"),
                            thr = NULL, min = NULL) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  type <- match.arg(type)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }

  if(type == "complete") {
    keep <- !apply(assay(se), 1, function(x) any(is.na(x)))
    filtered <- se[keep,]
  }
  if(type == "condition") {
    assertthat::assert_that(is.numeric(thr),
                            length(thr) == 1)
    max_repl <- max(colData(se)$replicate)
    if(thr < 0 | thr > max_repl) {
      stop("invalid filter threshold 'thr' applied",
           "\nRun filter() with a threshold ranging from 0 to ",
           max_repl)
    }

    filtered <- filter_missval(se, thr = thr)
  }
  if(type == "fraction") {
    assertthat::assert_that(is.numeric(min),
                            length(min) == 1)
    if(min < 0 | min > 1) {
      stop("invalid filter threshold 'min' applied",
           "\nRun filter() with a percent ranging from 0 to 1")
    }

    bin_data <- assay(se)
    idx <- is.na(assay(se))
    bin_data[!idx] <- 1
    bin_data[idx] <- 0

    # Filter se on the maximum allowed number of
    # missing values per condition (defined by thr)
    keep <- bin_data %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      gather(ID, value, -rowname) %>%
      group_by(rowname) %>%
      summarize(n = n(),
                valid = sum(value),
                frac = valid / n) %>%
      filter(frac >= min)
    filtered <- se[keep$rowname, ]
  }
  return(filtered)
}

#' Normalization using vsn
#'
#' \code{normalize_vsn} performs variance stabilizing transformation
#' using the \code{\link[vsn]{vsn-package}}.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()}.
#' @return A normalized SummarizedExperiment object.
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
#' # Filter and normalize
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' @export
normalize_vsn <- function(se) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # Variance stabilization transformation on assay data
  se_vsn <- se
  vsn.fit <- vsn::vsnMatrix(2 ^ assay(se_vsn))
  assay(se_vsn) <- vsn::predict(vsn.fit, 2 ^ assay(se_vsn))
  return(se_vsn)
}

#' Imputation by random draws from a manually defined distribution
#'
#' \code{manual_impute} imputes missing values in a proteomics dataset
#' by random draws from a manually defined distribution.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()}
#' and normalize the data using \code{\link{normalize_vsn}()}.
#' @param shift Numeric(1),
#' Sets the left-shift of the distribution (in standard deviations) from
#' the median of the original distribution.
#' @param scale Numeric(1),
#' Sets the width of the distribution relative to the
#' standard deviation of the original distribution.
#' @return An imputed SummarizedExperiment object.
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
#' # Filter and normalize
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#'
#' # Impute missing values manually
#' imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
#' @export
manual_impute <- function(se, scale = 0.3, shift = 1.8) {
  if(is.integer(scale)) scale <- is.numeric(scale)
  if(is.integer(shift)) shift <- is.numeric(shift)
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(scale),
                          length(scale) == 1,
                          is.numeric(shift),
                          length(shift) == 1)

  se_assay <- assay(se)

  # Show error if there are no missing values
  if(!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), "'",
         call. = FALSE)
  }

  # Get descriptive parameters of the current sample distributions
  stat <- se_assay %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(samples, value, -rowname) %>%
    filter(!is.na(value))  %>%
    group_by(samples) %>%
    summarise(mean = mean(value),
              median = median(value),
              sd = sd(value),
              n = n(),
              infin = nrow(se_assay) - n)
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
  for (a in seq_len(nrow(stat))) {
    assay(se)[is.na(assay(se)[, stat$samples[a]]), stat$samples[a]] <-
      rnorm(stat$infin[a],
            mean = stat$median[a] - shift * stat$sd[a],
            sd = stat$sd[a] * scale)
  }
  return(se)
}

#' Deprecated Function to coerce SummarizedExperiment to MSnSet object
#'
#' Use \code{\link[methods]{as}} instead.
#'
#' @param se SummarizedExperiment,
#' Object which will be turned into a MSnSet object.
#' @return A MSnSet object.
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
#' # Convert to MSnSet
#' data_msn <- as(se, "MSnSet")
#' # Convert back to SE
#' se_back <- as(data_msn, "SummarizedExperiment")
#' @export
se2msn <- function(se) {
  .Deprecated("as(se, \"MSnSet\")")
  as(se, "MSnSet")
}

#' Impute missing values
#'
#' \code{impute} imputes missing values in a proteomics dataset.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()}
#' and normalize the data using \code{\link{normalize_vsn}()}.
#' @param fun "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "man", "min", "zero", "mixed" or "nbavg",
#' Function used for data imputation based on \code{\link{manual_impute}}
#' and \code{\link[MSnbase:impute-methods]{impute}}.
#' @param ... Additional arguments for imputation functions as depicted in
#' \code{\link{manual_impute}} and \code{\link[MSnbase:impute-methods]{impute}}.
#' @return An imputed SummarizedExperiment object.
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
#' # Filter and normalize
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#'
#' # Impute missing values using different functions
#' imputed_MinProb <- impute(norm, fun = "MinProb", q = 0.05)
#' imputed_QRILC <- impute(norm, fun = "QRILC")
#'
#' imputed_knn <- impute(norm, fun = "knn", k = 10, rowmax = 0.9)
#' imputed_MLE <- impute(norm, fun = "MLE")
#'
#' imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
#' @export
impute <- function(se, fun = c("bpca", "knn", "QRILC", "MLE",
                               "MinDet", "MinProb", "man", "min", "zero",
                               "mixed", "nbavg"), ...) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(fun))

  # Show error if inputs do not contain required columns
  fun <- match.arg(fun)

  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }

  # Show error if there are no missing values
  if(!any(is.na(assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)), "'. ",
            "Returning the unchanged object.",
         call. = FALSE)
    return(se)
  }

  # Annotate whether or not there are missing values and how many
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))

  # if the "man" function is selected, use the manual impution method
  if(fun == "man") {
    se <- manual_impute(se, ...)
  }
  # else use the MSnSet::impute function
  else {
    MSnSet_data <- as(se, "MSnSet")
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
    assay(se) <- MSnbase::exprs(MSnSet_imputed)
  }
  return(se)
}

#' Differential enrichment test
#'
#' \code{test_diff} performs a differential enrichment test based on
#' protein-wise linear models and empirical Bayes
#' statistics using \pkg{limma}. False Discovery Rates are estimated
#' using \pkg{fdrtool}.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()},
#' normalize the data using \code{\link{normalize_vsn}()} and
#' impute remaining missing values using \code{\link{impute}()}.
#' @param type "control", "all" or "manual",
#' The type of contrasts that will be tested.
#' This can be all possible pairwise comparisons ("all"),
#' limited to the comparisons versus the control ("control"), or
#' manually defined contrasts ("manual").
#' @param control Character(1),
#' The condition to which contrasts are generated if type = "control"
#' (a control condition would be most appropriate).
#' @param test Character,
#' The contrasts that will be tested if type = "manual".
#' These should be formatted as "SampleA_vs_SampleB" or
#' c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @param design_formula Formula,
#' Used to create the design matrix.
#' @return A SummarizedExperiment object
#' containing fdr estimates of differential expression.
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
#' diff <- test_diff(imputed, "manual",
#'     test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))
#'
#' # Test for differentially expressed proteins with a custom design formula
#' diff <- test_diff(imputed, "control", "Ctrl",
#'     design_formula = formula(~ 0 + condition + replicate))
#' @export
test_diff <- function(se, type = c("control", "all", "manual"),
                      control = NULL, test = NULL,
                      design_formula = formula(~ 0 + condition)) {

  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")

  # Show error if inputs do not contain required columns
  type <- match.arg(type)

  col_data <- colData(se)
  raw <- assay(se)

  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if(any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }

  if(!is.null(control)) {
    # Show error if control input is not valid
    assertthat::assert_that(is.character(control),
                            length(control) == 1)
    if(!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"), "'",
           call. = FALSE)
    }
  }

  # variables in formula
  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]

  # Throw error if variables are not col_data columns
  if(any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if(variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }

  # Obtain variable factors
  for(var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }

  # Make an appropriate design matrix
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))

  # Generate contrasts to be tested
  # Either make all possible combinations ("all"),
  # only the contrasts versus the control sample ("control") or
  # use manual contrasts
  conditions <- as.character(unique(condition))
  if(type == "all") {
    # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")

    if(!is.null(control)) {
      # Make sure that contrast containing
      # the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if(length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }

  }
  if(type == "control") {
    # Throw error if no control argument is present
    if(is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")

    # Make contrasts
    cntrst <- paste(conditions[!conditions %in% control],
                    control,
                    sep = " - ")
  }
  if(type == "manual") {
    # Throw error if no test argument is present
    if(is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))

    if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
           "'.", call. = FALSE)
    }

    cntrst <- gsub("_vs_", " - ", test)

  }
  # Print tested contrasts
  message("Tested contrasts: ",
          paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))

  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- contrasts.fit(fit, made_contrasts)

  if(any(is.na(raw))) {
    for(i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
    }
  }

  eB_fit <- eBayes(contrast_fit)

  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", coef = comp,
                    number = Inf, confint = TRUE)
    res <- res[!is.na(res$t),]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }

  # Retrieve the differential expression test results
  limma_res <- map_df(cntrst, retrieve_fun)

  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    select(rowname, logFC, CI.L, CI.R, P.Value, qval, comparison) %>%
    mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    gather(variable, value, -c(rowname,comparison)) %>%
    mutate(variable = recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
    unite(temp, comparison, variable) %>%
    spread(temp, value)
  rowData(se) <- merge(rowData(se, use.names = FALSE), table,
    by.x = "name", by.y = "rowname", all.x = TRUE, sort=FALSE)
  return(se)
}

#' Mark significant proteins
#'
#' \code{add_rejections} marks significant proteins based on defined cutoffs.
#'
#' @param diff SummarizedExperiment,
#' Proteomics dataset on which differential enrichment analysis
#' has been performed (output from \code{\link{test_diff}()}).
#' @param alpha Numeric(1),
#' Sets the threshold for the adjusted P value.
#' @param lfc Numeric(1),
#' Sets the threshold for the log2 fold change.
#' @return A SummarizedExperiment object
#' annotated with logical columns indicating significant proteins.
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
#' @export
add_rejections <- function(diff, alpha = 0.05, lfc = 1) {
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1)

  row_data <- rowData(diff, use.names = FALSE) %>%
    as.data.frame()
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(diff)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(diff)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }

  # get all columns with adjusted p-values and log2 fold changes
  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))

  # Mark differential expressed proteins by
  # applying alpha and log2FC parameters per protein
  if(length(cols_p) == 1) {
    rowData(diff)$significant <-
      row_data[, cols_p] <= alpha & abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <-
      rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <-
      gsub("p.adj", "significant", colnames(row_data)[cols_p])
  }
  if(length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df,
      significant = apply(sign_df, 1, function(x) any(x)))
    colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))

    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    rowData(diff) <- merge(rowData(diff, use.names = FALSE), sign_df,
                           by = "name")
  }
  return(diff)
}

#' Generate a results table
#'
#' \code{get_results} generates a results table from a proteomics dataset
#' on which differential enrichment analysis was performed.
#'
#' @param dep SummarizedExperiment,
#' Data object for which differentially enriched proteins are annotated
#' (output from \code{\link{test_diff}()} and \code{\link{add_rejections}()}).
#' @return A data.frame object
#' containing all results variables from the performed analysis.
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
#' # Get results
#' results <- get_results(dep)
#' colnames(results)
#'
#' significant_proteins <- results[results$significant,]
#' nrow(significant_proteins)
#' head(significant_proteins)
#' @export
get_results <- function(dep) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))

  row_data <- rowData(dep, use.names = FALSE)
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }

  # Obtain average protein-centered enrichment values per condition
  row_data$mean <- rowMeans(assay(dep), na.rm = TRUE)
  centered <- assay(dep) - row_data$mean
  centered <- data.frame(centered) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(dep)), by = "ID")
  centered <- group_by(centered, rowname, condition) %>%
    summarize(val = mean(val, na.rm = TRUE)) %>%
    mutate(val = signif(val, digits = 3)) %>%
    spread(condition, val)
  colnames(centered)[2:ncol(centered)] <-
    paste(colnames(centered)[2:ncol(centered)], "_centered", sep = "")

  # Obtain average enrichments of conditions versus the control condition
  ratio <- as.data.frame(row_data) %>%
    column_to_rownames("name") %>%
    select(ends_with("diff")) %>%
    signif(., digits = 3) %>%
    rownames_to_column()
  colnames(ratio)[2:ncol(ratio)] <-
    gsub("_diff", "_ratio", colnames(ratio)[2:ncol(ratio)])
  df <- left_join(ratio, centered, by = "rowname")

  # Select the adjusted p-values and significance columns
  pval <- as.data.frame(row_data) %>%
    column_to_rownames("name") %>%
    select(ends_with("p.val"),
      ends_with("p.adj"),
      ends_with("significant")) %>%
    rownames_to_column()
  pval[, grep("p.adj", colnames(pval))] <-
    pval[, grep("p.adj", colnames(pval))] %>%
    signif(digits = 3)

  # Join into a results table
  ids <- as.data.frame(row_data) %>% select(name, ID)
  table <- left_join(ids, pval, by = c("name" = "rowname"))
  table <- left_join(table, df, by = c("name" = "rowname")) %>%
    arrange(desc(significant))
  return(table)
}

#' Generate a wide data.frame from a SummarizedExperiment
#'
#' \code{get_df_wide} generate a wide data.frame from a SummarizedExperiment.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}).
#' @return A data.frame object
#' containing all data in a wide format,
#' where each row represents a protein.
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

#' Generate a long data.frame from a SummarizedExperiment
#'
#' \code{get_df_long} generate a wide data.frame from a SummarizedExperiment.
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}).
#' @return A data.frame object
#' containing all data in a wide format,
#' where each row represents a single measurement.
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
#' # Get a long data.frame
#' long <- get_df_long(dep)
#' colnames(long)
#' @export
get_df_long <- function(se) {
  # Show error if inputs are not the required classes
  assert_that(inherits(se, "SummarizedExperiment"))

  # Show error if inputs do not contain required columns
  if (!"name" %in% colnames(rowData(se, use.names = FALSE))) {
    stop("'name' column is not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }

# Extract column data
col_data <- colData(se) %>%
  data.frame() %>%
  rownames_to_column() %>%
  select(-ID)
# Extract row data
row_data <- rowData(se, use.names = FALSE) %>%
  data.frame()
# Extract assay data
assay_data <- assay(se) %>%
  data.frame() %>%
  rownames_to_column()
colnames(assay_data)[1] <- "name"

# Transform assay_data in long format
long_assay <- assay_data %>%
  gather("rowname", "intensity", -name)

# Merge row and assay data into a wide data.frame
long <- long_assay %>%
  full_join(col_data, ., by = "rowname") %>%
  select(-rowname) %>%
  full_join(., row_data, by = "name")

return(long)
}
