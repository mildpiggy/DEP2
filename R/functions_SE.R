#' Transform NA to 0 in assay of a SummarizedExperiment
#'
#' @param se A SummarizedExperiment object which contain NAs in assay.
#'
#' @return
#' A SummarizedExperiment object
#' @export
#'
NAiszero <- function(se){
  # print(class(assay(se)[1,1]))
  assay(se)[is.na(assay(se))] = as.integer(0)
  # print(class(assay(se)[1,1]))
  return(se)
}



#' Filter SummarizedExperiment on missing values or formula
#'
#' \code{filter_se} filters a SummarizedExperiment object based on missing values in assay
#' and provides formula.
#' The dataset is filtered for proteins that have a maximum of
#' 'thr' missing values in at least one condition; or a maximum of 'missnum' in all samples;
#' or a filter_formula
#'
#' @param se SummarizedExperiment,
#' Proteomics data (output from \code{\link{make_se}()} or
#' \code{\link{make_se_parse}()}).
#' @param thr Integer(1),
#' Sets the threshold for the allowed max number of missing values
#' in at least one condition. At least in one condition, missing values number < or = 'thr'.
#' @param missnum Integer(1),
#' Sets the threshold for the allowed max number among all samples,
#' missing values number < or = 'missnum' in each row.
#' @param fraction A numeric from 0 to 1, threshold of missing occupancy of each row
#' @param filter_formula Formula or character(1).
#' A filter expression, written as a formula.
#' @param rowsum_threshold A numeric, row sum of intensity(abundance) must larger than it.
#' @param return_keeprows Logical(1), default FALSE, return the trimmed object.
#' If TURE return the rows number that pass through filter.
#' @return A filtered SummarizedExperiment object,
#' or a vector of keep rows if \code{return_keeprows} is TRUE
#' @examples
#' # Load example
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")
#'
#' # Filter
#' stringent_filter <- filter_se(se, thr = 0, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' less_stringent_filter <- filter_se(se, thr = 1, missnum = 2, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' keep_rows <- filter_se(se, thr = 0, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+", return_keeprows = TRUE)
#' @export
filter_se <- function(se,
                      thr = NULL,
                      missnum = NULL,
                      fraction = NULL,
                      filter_formula = NULL,
                      rowsum_threshold = NULL,
                      return_keeprows = FALSE){
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.null(thr)|is.numeric(thr)|is.integer(thr), is.null(missnum)|is.numeric(missnum)|is.integer(missnum),
                          is.null(filter_formula)|is.character(filter_formula)|purrr::is_formula(filter_formula),
                          is.null(fraction)|is.numeric(fraction), is.null(fraction)|length(fraction) == 1)
  origin_se = se

  if(inherits(se,"DESeqDataSet")){
    se <- QFeatures::zeroIsNA(se) ## DESeqDataSet missing values is store as zero
  }
  if(!is.null(missnum)){
    message("filter base on missing number <= ",missnum,"\n")
    rowData(se)$nNonZero <- rowSums(!is.na(assay(se)))
    keeprow = (rowData(se)$nNonZero >= (ncol(assay(se)) - missnum))
    se <- se[keeprow,]
  }
  if(!is.null(thr)){
    message("filter base on missing number is <= ",thr," in at least one condition.\n")
    se <- filter_missval(se, thr)
  }

  if(!is.null(fraction)){
    message("filter base on missing number fraction < ",fraction," in each column\n")
    if(fraction < 0 | fraction > 1) {
      warning("invalid filter threshold 'fraction' applied",
              "Run filter_se() with a percent ranging from 0 to 1.\n",
              "Skipped the filter by 'fraction'")
    }else{
      se <- filter_fraction(se, fraction)
    }
  }

  if(!is.null(filter_formula)){
    message("filter base on giving formula \n")
    se <- filterSummarizedExperimentWithFormula(object = se, filter = filter_formula)
  }

  if(!is.null(rowsum_threshold)){
    message("filter base on counts sum >=",rowsum_threshold,"\n")
    se = se[rowSums(assay(se),na.rm = T) >= rowsum_threshold, ]
  }
  if(return_keeprows)
    return(which(rownames(origin_se) %in% rownames(se)))


  if(inherits(se,"DESeqDataSet")){
    se <- NAiszero(se) ## DESeqDataSet missing values is store as zero
  }

  return(se)
}


# filter on the missing number threashold in one condition
filter_missval <- function(se, thr = 0) {
  # Show error if inputs are not the required classes
  if(is.integer(thr)) thr <- as.numeric(thr)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(thr),
                          length(thr) == 1)

  if(any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  max_repl <- max(colData(se)$replicate %>% as.numeric())
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
  suppressMessages(keep <- bin_data %>%
                     data.frame() %>%
                     rownames_to_column() %>%
                     gather(ID, value, -rowname) %>%
                     left_join(., data.frame(colData(se)), by = "ID") %>%
                     group_by(rowname, condition) %>%
                     summarize(miss_val = n() - sum(value)) %>%
                     filter(miss_val <= thr) %>%
                     spread(condition, miss_val))

  se_fltrd <- se[keep$rowname, ]
  return(se_fltrd)
}

filter_fraction <- function(se, fraction){
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
    filter(frac >= fraction)
  filtered <- se[keep$rowname, ]
}

##' @importFrom lazyeval f_eval
filterSummarizedExperimentWithFormula <- function(object, filter, na.rm = TRUE,...) {
  x <- rowData(object)
  sel <- tryCatch(lazyeval::f_eval(filter, data = as.list(x)),
                  error = function(e) rep(FALSE, nrow(x)))
  sel[is.na(sel)] <- !na.rm
  object[sel,  ]
}
