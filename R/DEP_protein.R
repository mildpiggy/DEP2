

#' Clean characters in expression columns
#'
#' If the class of expression contain characters(e.g. 'NA.', 'filtered'), run Clean characters to replace character
#' by NA, and transfor expression value to numeric
#'
#' @param express_assay data.frame or matrix, the expression assay.
#'
#' @return
#' a numeric matrix
#' @export
#'
#' @examples
#' expression_assay = data.frame(A = c("1165799.5","51419028","1.9e6","2213111.25"),
#'                               B = c("1455111","Filtered","2021020.375","NA"),
#'                               c = c(261344, 49840292, 2388750, 1753475.625))
#' str(expression_assay)
#' (clean_character(expression_assay))
#'
clean_character <- function(express_assay){
  # assertthat::assert_that(is.data.frame(express_assay)|is.matrix(express_assay))
  express_assay <- as.data.frame(express_assay) %>%
    apply(., 2, function(x){suppressWarnings(as.numeric(x))}) %>% as.data.frame()
}



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
#' @importFrom dplyr mutate
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
  if(is.null(names) & is.null(ids)){
    stop("Either names or ids must be provide.")
  }
  if(is.null(names)){
    warning("Names do not input and use ids instead.")
    names = ids
  }else if(is.null(ids)){
    warning("Ids do not input and use names instead.")
    ids = names
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
    dplyr::mutate(name = gsub(paste0(delim, ".*"), "", get(names)),
                  ID = gsub(paste0(delim, ".*"), "", get(ids)),
                  name = make.unique(ifelse(name == "" |
                                              is.na(name), ID, name)))
  return(proteins_unique)
}

# same to DEP
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
#' @param log2transform Logical(1), whether log2 transform the assay, default TRUE.
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
make_se <- function (proteins_unique, columns, expdesign, log2transform = TRUE)
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.numeric(columns)|is.integer(columns)| is.character(columns),
                          is.data.frame(expdesign))
  if(is.numeric(columns)) columns = as.integer(columns)
  if(is.character(columns)){
    if(!all(columns %in% colnames(proteins_unique)))
      stop("columns should be the columns in ", deparse(substitute(proteins_unique)), "but ", columns[!columns %in% colnames(proteins_unique)],"do not exist.")
    columns = which(colnames(proteins_unique) %in% columns)
  }else if(is.integer(columns)){
    assert_that(all(columns %in% 1:nrow(proteins_unique)))
  }
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design", call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se() with the appropriate columns as argument",
         "Or run clean_character() to filter character values(like 'NA'), and transfer numeric in these columns",
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if (tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA  ## zero to NA

  if(log2transform){
    raw <- log2(raw)
  }
  expdesign <- mutate(expdesign, label = make.names(label), condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))) )
  if(any(is.na(matched))
     # && sum(matched,na.rm = T) != length(ecols)
  ){
    matched <- match(make.names(expdesign$label),
                     make.names(colnames(raw)) )
  }
  if (any(is.na(matched))
      # && sum(matched,na.rm = T) != length(ecols)
  ) {
    stop("Labels of the experimental design do not match ",
         "with column names in 'proteins_unique'", "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]

  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  se <- SummarizedExperiment(assays = as.matrix(raw), colData = expdesign,
                             rowData = row_data)
  return(se)
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
#' @param remove_prefix Logical(1),
#' whether remove the prefix of expression columns.
#' @param log2transform Logical(1), whether log2 transform the assay, default TRUE.
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
make_se_parse <- function (proteins_unique, columns, mode = c("char", "delim"),
                          chars = 1, sep = "_", remove_prefix = T, remove_suffix = F, log2transform = T)
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.numeric(columns)|is.integer(columns)| is.character(columns),
                          is.character(mode), is.numeric(chars), length(chars) ==
                            1, is.character(sep), length(sep) == 1)
  if(is.numeric(columns)) columns = as.integer(columns)
  if(is.character(columns)){
    if(!all(columns %in% colnames(proteins_unique)))
      stop("columns should be the columns in ", deparse(substitute(proteins_unique)), "but ", columns[!columns %in% colnames(proteins_unique)],"do not exist.")
    columns = which(colnames(proteins_unique) %in% columns)
  }else if(is.integer(columns)){
    assert_that(all(columns %in% 1:nrow(proteins_unique)))
  }
  mode <- match.arg(mode)
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("Specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         "Or run clean_character() to filter character values(like 'NA'), and transfer numeric in these columns",
         call. = FALSE)
  }

  if (tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  if(log2transform){
    raw <- log2(raw)
  }
  if(remove_prefix){
    colnames(raw) <- delete_prefix(colnames(raw)) %>% make.names()
  }
  if(remove_suffix){
    colnames(raw) <- delete_suffix(colnames(raw)) %>% make.names()
  }
  row_data <- proteins_unique[, -columns]

  rownames(row_data) <- row_data$name
  if (mode == "char") {
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      mutate(condition = substr(label, 1, nchar(label) -
                                  chars), replicate = substr(label, nchar(label) +
                                                               1 - chars, nchar(label))) %>% unite(ID, condition,
                                                                                                   replicate, remove = FALSE)
  }
  if (mode == "delim") {
    # colnames(raw) = gsub(get_suffix(colnames(raw)),"", colnames(raw))
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      separate(label, c("condition", "replicate"), sep = sep,
               remove = FALSE, extra = "merge") %>% unite(ID,
                                                          condition, replicate, remove = FALSE)
  }
  rownames(col_data) <- col_data$ID
  colnames(raw)[match(col_data$label, colnames(raw))] <- col_data$ID
  raw <- raw[, !is.na(colnames(raw))]
  se <- SummarizedExperiment(assays = as.matrix(raw), colData = col_data,
                             rowData = row_data)
  return(se)
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
normalize_vsn <- function (se)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_vsn <- se
  vsn.fit <- vsn::vsnMatrix(2^assay(se_vsn))
  assay(se_vsn) <- vsn::predict(vsn.fit, 2^assay(se_vsn))
  return(se_vsn)
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
#' "MinProb", "man", "min", "zero", "mixed", "nbavg", "GSimp" or "RF",
#' Function used for data imputation based on \code{\link{manual_impute}}("man")
#' , \code{\link[MSnbase]{impute,MSnSet-method}} and \code{\link[missForest]{missForest}}("RF").
#' @param ... Additional arguments for imputation functions as depicted in
#' \code{\link{manual_impute}}, \code{\link[missForest]{missForest}} and [MsCoreUtils::impute_matrix()].
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
#' @importFrom MSnbase impute exprs
#' @importFrom missForest missForest
impute <- function (se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet",
                                "MinProb", "man", "min", "zero", "mixed", "nbavg","RF", "GSimp"), ...)
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(fun))
  fun <- match.arg(fun)
  # if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
  #   stop("'name' and/or 'ID' columns are not present in '",
  #        deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns",
  #        call. = FALSE)
  # }
  if (!any(is.na(assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)),
            "'. ", "Returning the unchanged object.", call. = FALSE)
    return(se)
  }
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))

  if (fun == "man") {
    se <- manual_impute(se, ...)
  }else if(fun == "RF"){
    args <- list(...)
    if(! "ntree" %in% names(args))
      args[["ntree"]] = 60
    if(! "parallelize" %in% names(args) )
      args[["parallelize"]] = "variables"
    args$xmis = assay(se)
    doParallel::registerDoParallel(cores=2)
    assay(se) <- do.call(missForest::missForest, args)$ximp
    doParallel::stopImplicitCluster()
  }else if(fun == "GSimp"){
    assay(se) <- GS_imp_wrapper(assay(se), ...)
  # }else if(fun == "GMSimpute"){
  #   impassay <- GMSimpute::GMS.Lasso(assay(se),log.scale=F,TS.Lasso=TRUE)
  #   cat("finish")
  #   se = se[match(rownames(impassay), rownames(se))]
  #   assay(se) <- impassay
  }else{
    MSnSet_data <- as(se, "MSnSet") # transfer to MSnSet object
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun,
                                      ...)
    assay(se) <- MSnbase::exprs(MSnSet_imputed)
  }
  return(se)
}

## This function un-export
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


# setGeneric("add_rejections", function(diff, alpha = 0.05, lfc = 1,thresholdmethod="intersect",curvature=1,x0_fold = 2,...) {
#   standardGeneric("add_rejections")
# })

#' Differential enrichment/expression test
#'
#' \code{test_diff} performs a differential enrichment/expression test based on
#' protein/peptide-wise linear models and empirical Bayes
#' statistics using \pkg{limma}. False Discovery Rates are estimated
#' using \code{"Strimmer's qvalue"}, \code{"Benjamini-Hochberg fdr"} or \code{"Storey's qvalue"}.
#'
#' @param se SummarizedExperiment,
#' Proteomics data from quantity table (output from \code{\link{make_se}()},
#' \code{\link{make_se_parse}()}), or from re-aggregate peptide quantity from \code{\link{aggregateFeatures}()} ).
#' For proteomics data directly from quantity table, it is adviced to first remove
#' low quanlity and undesired rows by \code{\link{filter_se}()}, and missing values are imputed by \code{\link{impute}()}.
#' Normalization is recommended in both workflow.
#'
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
#' These should be formatted as "conditionA_vs_conditionB" or
#' c("conditionA_vs_conditionC", "conditionB_vs_conditionC").
#' @param design_formula Formula,
#' Used to create the design matrix.
#' @param fdr.type Character(1), the method to control false discovery rate,
#' one of "Strimmer's qvalue(t)","Strimmer's qvalue(p)","BH","Storey's qvalue".
#' "Strimmer's qvalue" calculate fdr vis \code{\link[fdrtool]{fdrtool}()} using t-statistic or p values from limma.
#' "BH" calculate fdr vis \code{\link[stats]{p.adjust}()} using "BH" method .
#' "Storey's qvalue"  calculate fdr vis \code{\link{stats::qvalue}()}.
#' @param contrast_upon C the contrast to bulid upon which column in experimentdesign
#' @param filter_ihw logical(1),Whether to use ihw to perform independent filtering and p-value adjustment, default FALSE
#' @param independentFiltering logical(1), whether independent filtering should be applied automatically, default TRUE
#' @param lfcshark logical(1),whether to use lfcshark, default FALSE
#'
#' @return A SummarizedExperiment object or a DEGdata
#' containing diff, test result from limma and fdr estimates of differential expression.
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
#' @importFrom fdrtool fdrtool
#' @export
test_diff <- function(se, type = c("all", "control", "manual"),
                      control = NULL, test = NULL,
                      design_formula = formula(~ 0 + condition),
                      fdr.type = c("Strimmer's qvalue(t)","Strimmer's qvalue(p)","BH","Storey's qvalue")) {

  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.null(control) | is.character(control),
                          is.null(test) | is.character(test),
                          # is.character(contrast_upon),
                          is.character(type),
                          class(design_formula) == "formula",
                          is.character(fdr.type))
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  fdr.type <- match.arg(fdr.type)

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

  # Show error if inputs do not contain required columns
  # fdr.type <- match.arg(fdr.type)
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

  eB_fit <- eBayes(contrast_fit,trend = FALSE)

  # function to retrieve the results of
  # the differential expression test using 'fdrtool',('BH','qvalue')
  retrieve_fun <- function(comp, fit = eB_fit, fdr.type){
    res <- topTable(fit, sort.by = "t", coef = comp,
                    number = Inf, confint = TRUE)
    res <- res[!is.na(res$t),]

    if(fdr.type == "Strimmer's qvalue(t)"){ ##
      fdr_res <- fdrtool(res$t,plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
    }
    if(fdr.type == "Strimmer's qvalue(p)"){
      fdr_res <- fdrtool(res$P.Value,statistic = "pvalue",plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
    }
    if(fdr.type == "BH"){
      padj <- p.adjust(res$P.Value,method = "BH")
      res$qval = padj
    }
    if(fdr.type == "Storey's qvalue"){
      qval_res = qvalue::qvalue(res$P.Value)
      res$qval = qval_res$qvalues
    }

    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }

  message(fdr.type)
  # Retrieve the differential expression test results
  limma_res <- purrr::map_df(cntrst, retrieve_fun,fdr.type = fdr.type)

  table <- limma_res %>%
    select(rowname, logFC, CI.L, CI.R, t,P.Value, qval, comparison) %>%
    mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    gather(variable, value, -c(rowname,comparison)) %>%
    mutate(variable = recode(variable, logFC = "diff", t="t.stastic", P.Value = "p.val", qval = "p.adj")) %>%
    unite(temp, comparison, variable) %>%
    spread(temp, value)

  rowData(se) <- merge(rowData(se, use.names = FALSE), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE, sort=FALSE)
  return(se)
}


# diff <- deg_save
# add_rejections <- function(object, alpha = 0.05, lfc = 1,thresholdmethod="intersect",curvature=1,x0_fold = 2,...){
#   UseMethod('add_rejections')
# }



setGeneric("add_rejections", function(diff, alpha = 0.05, lfc = 1,thresholdmethod="intersect",curvature=1,x0_fold = 2,...) {
  standardGeneric("add_rejections")
})


#' Mark significant proteins
#'
#' \code{add_rejections} marks significant proteins based on defined cutoffs.
#'
#' @param diff SummarizedExperiment,
#' Proteomics dataset on which differential analysis
#' has been performed (output from \code{\link{test_diff}()}) or
#' RNA-seq dataset (output from \code{\link{}get_res}()})
#' @param alpha Numeric(1),
#' Sets the threshold for the adjusted P value.
#' @param lfc Numeric(1),
#' Sets the threshold for the log2 fold change.
#' @param thresholdmethod Character(1),
#' The thresholdmethod to decide significant.should be one of "intersect" ,"curve". If thresholdmethod = "curve", cutoff lines is the curves with function y = c/(x-x0). c is the curvature, and x0 is the minimum fold change
#' @param curvature Numeric(1),
#' Sets the curvature for the curve cutoff lines
#' @param x0_fold Numeric(1),
#' decide the x0 ratio to the standard deviations of L2FC. The x0 usually is set to 1(medium confidence) or 2(high confidence) standard deviations.
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
#' @name add_rejections
#' @rdname DEP2_add_rejections
NULL

##' @exportMethod add_rejections
##' @rdname DEP2_add_rejections
setMethod("add_rejections",
          "SummarizedExperiment",
          function(diff, ...){
            add_rejections.SummarizedExperiment(diff = diff,
                                                alpha = alpha, lfc = lfc, thresholdmethod = thresholdmethod, curvature = curvature, x0_fold = x0_fold
                                                # ...
                                                )
          }
)
# add_rejections(my_dep_save)


add_rejections.SummarizedExperiment <- function (diff, alpha = 0.05, lfc = 1,thresholdmethod="intersect",curvature=1,x0_fold = 2)
{
  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(diff)), "'\nRun make_unique() and make_se()/make_pe() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }

  if( length(thresholdmethod)!=1 || (!thresholdmethod %in% c("intersect","curve")) ){
    stop("thresholdmethod should be one of intersect/curve")
  }

  if( length(grep("vs(.)*_significant$|^significant$", colnames(row_data))) > 0 ){
    exit_sig_cols <- grep("_significant|significant", colnames(row_data))
    warning("Significant columns already exited in row_data.\n   Remove/replace values in under columns",
            paste(colnames(row_data)[exit_sig_cols],collapse ="; ") )
    rowData(diff) = row_data = row_data[,-exit_sig_cols]
  }

  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))

  if(thresholdmethod == "intersect"){
    if (is.integer(alpha))
      alpha <- as.numeric(alpha)
    if (is.integer(lfc))
      lfc <- as.numeric(lfc)
    assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                            is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                            length(lfc) == 1)
    if (length(cols_p) == 1) {
      rowData(diff)$significant <- row_data[, cols_p] <= alpha &
        abs(row_data[, cols_diff]) >= lfc
      rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
      colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub("p.adj",
                                                                              "significant", colnames(row_data)[cols_p])
    }
    if (length(cols_p) > 1) {
      p_reject <- row_data[, cols_p] <= alpha
      p_reject[is.na(p_reject)] <- FALSE
      diff_reject <- abs(row_data[, cols_diff]) >= lfc
      diff_reject[is.na(diff_reject)] <- FALSE
      sign_df <- p_reject & diff_reject
      sign_df <- cbind(sign_df, significant = apply(sign_df,
                                                    1, function(x) any(x)))
      colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
      sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
      rowData(diff) <- left_join(rowData(diff, use.names = FALSE) %>%as.data.frame(),
                                 sign_df, by = "name")

    }
    return(diff)

  }else if(thresholdmethod == "curve"){
    if (is.integer(curvature))
      curvature <- as.numeric(curvature)
    if (is.integer(x0_fold))
      x0_fold <- as.numeric(x0_fold)
    assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                            is.numeric(curvature), length(curvature) == 1, is.numeric(x0_fold),
                            length(x0_fold) == 1)

    if (length(cols_p) == 1) {
      x = row_data[, cols_diff]
      x = fun.outlier(x) %>% stats::na.omit()
      fit <- fitnormal(x)
      σ = sqrt(fit$theta[2])
      x0 = x0_fold*σ
      polar = ifelse((rowData(diff)[,cols_diff] > 0),1,-1)
      rowData(diff)$significant <- ( -log10(rowData(diff)[,cols_p]) > curvature/abs(rowData(diff)[,cols_diff] - polar*x0) &
                                       ifelse((row_data[,cols_diff] >= 0),row_data[,cols_diff]>x0,row_data[,cols_diff] < -x0) )
      rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
      colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub("p.adj",
                                                                              "significant", colnames(row_data)[cols_p])
    }
    if(length(cols_p) > 1){
      sign_df <- sapply(1:length(cols_p), function(i) {
        x = row_data[, cols_diff[i]]
        x = fun.outlier(x) %>% na.omit()
        fit <- fitnormal(x)
        σ = sqrt(fit$theta[2])
        x0 = x0_fold*σ
        polar = ifelse((rowData(diff)[,cols_diff[i]] > 0),1,-1)
        contrast_significant <- ( -log10(rowData(diff)[,cols_p[i]]) > curvature/abs(rowData(diff)[,cols_diff[i]] - polar*x0) &
                                    ifelse((row_data[,cols_diff[i]] >= 0),row_data[,cols_diff[i]]>x0,row_data[,cols_diff[i]] < -x0) )
      })
      colnames(sign_df) =  gsub("p.adj", "significant", colnames(row_data)[cols_p])
      sign_df <- cbind(sign_df, significant = apply(sign_df,
                                                    1, function(x) any(x)))
      sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
      rowData(diff) <- left_join(rowData(diff, use.names = FALSE) %>% as.data.frame(),
                                 sign_df, by = "name")
    }
    return(diff)
  }else stop("thresholdmethod should be one of 'intersect'/'curve'",
             call. = FALSE)

}



## rm outlier
fun.outlier <- function(x,time.iqr=1.5) {
  outlier.low <- quantile(x,probs=c(0.25))-IQR(x)*time.iqr
  outlier.high <- quantile(x,probs=c(0.75))+IQR(x)*time.iqr
  x[which(x>outlier.high | x<outlier.low)]<-NA
  x
}

## fit Gaussian distribution for numeric vector x
fitnormal <- function (x, exact = TRUE) {
  if (exact) {
    ################################################
    ## Exact inference based on likelihood theory ##
    ################################################
    ## minimum negative log-likelihood (maximum log-likelihood) estimator of `mu` and `phi = sigma ^ 2`
    n <- length(x)
    mu <- sum(x) / n
    phi <- crossprod(x - mu)[1L] / n  # (a bised estimator, though)
    ## inverse of Fisher information matrix evaluated at MLE
    invI <- matrix(c(phi, 0, 0, phi * phi), 2L,
                   dimnames = list(c("mu", "sigma2"), c("mu", "sigma2")))
    ## log-likelihood at MLE
    loglik <- -(n / 2) * (log(2 * pi * phi) + 1)
    ## return
    return(list(theta = c(mu = mu, sigma2 = phi), vcov = invI, loglik = loglik, n = n))
  }
  else {
    ##################################################################
    ## Numerical optimization by minimizing negative log-likelihood ##
    ##################################################################
    ## negative log-likelihood function
    ## define `theta = c(mu, phi)` in order to use `optim`
    nllik <- function (theta, x) {
      (length(x) / 2) * log(2 * pi * theta[2]) + crossprod(x - theta[1])[1] / (2 * theta[2])
    }
    ## gradient function (remember to flip the sign when using partial derivative result of log-likelihood)
    ## define `theta = c(mu, phi)` in order to use `optim`
    gradient <- function (theta, x) {
      pl2pmu <- -sum(x - theta[1]) / theta[2]
      pl2pphi <- -crossprod(x - theta[1])[1] / (2 * theta[2] ^ 2) + length(x) / (2 * theta[2])
      c(pl2pmu, pl2pphi)
    }
    ## ask `optim` to return Hessian matrix by `hessian = TRUE`
    ## use "..." part to pass `x` as additional / further argument to "fn" and "gn"
    ## note, we want `phi` as positive so box constraint is used, with "L-BFGS-B" method chosen
    init <- c(sample(x, 1), sample(abs(x) + 0.1, 1))  ## arbitrary valid starting values
    z <- optim(par = init, fn = nllik, gr = gradient, x = x, lower = c(-Inf, 0), method = "L-BFGS-B", hessian = TRUE)
    ## post processing ##
    theta <- z$par
    loglik <- -z$value  ## flip the sign to get log-likelihood
    n <- length(x)
    ## Fisher information matrix (don't flip the sign as this is the Hessian for negative log-likelihood)
    I <- z$hessian / n  ## remember to take average to get mean
    invI <- solve(I, diag(2L))  ## numerical inverse
    dimnames(invI) <- list(c("mu", "sigma2"), c("mu", "sigma2"))
    ## return
    return(list(theta = theta, vcov = invI, loglik = loglik, n = n))
  }
}





#' Set order of experiment design for SummarizedExperiment
#'
#' \code{Order_cols} can re-set the order of a SummarizedExperiment object according  colData()
#'
#' @param object a SummarizedExperiment object
#' @param order  Character vector, the order to set.
#' @param order_column Set order on which column of experiment design, e.g. "condition", "label"
#'
#' @return
#' A object sorted by giving order
#' @export
#'
#' @examples
Order_cols <- function(object, order, order_column = "condition"){
  assertthat::assert_that(is.character(order), is.character(order_column) && length(order_column) == 1,
                          inherits(object, "SummarizedExperiment"))
  cd <- colData(object)
  if(nrow(cd) == 0 | ncol(cd) == 0){
    stop("object should contain a colData")
  }
  if(!order_column %in% colnames(cd)){
    stop("order_column should exist in the colData of object")
  }
  if(!(all(order %in% cd[,order_column]) && all(cd[,order_column] %in% order))){
    stop("order should include all factors in '", order_column, "' column: ",unique(paste0(cd[,order_column],seq =" ")) )
  }
  colData(object)[,order_column] = factor(colData(object)[,order_column], levels = unique(order))
  object = object[,order(colData(object)[,order_column])]
  return(object)
}








