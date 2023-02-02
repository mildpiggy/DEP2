

#' Data.frame to QFeatures object conversion using an experimental design
#'
#' \code{make_pe} creates a QFeatures object
#' based on two data.frames: the peptide table and experimental design.
#'
#' @param Peptide
#' Data.frame,
#' Peptide table with its belonging proteingroup.
#' @param columns A numeric indicating the indices of the columns
#' to be used as expression values. Can also be a character indicating the names of the columns.
#' @param expdesign A data.frame.
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' @param fnames An optional character(1) or numeric(1) indicating the column identifier.
#' @param assay_name
#' An character(1) to name assay in the QFeatures object.
#' @param log2transform Logical(1), whether log2 transform the assay, default TRUE.
#' @return
#' A QFeatures object with a peptide assay.
#' @export
#'
#' @examples
#' # Load a peptides
#' data(Silicosis_peptide)
#' # Exctract Intensity columns
#' (ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = T))
#' (expDesign = Silicosis_ExpDesign)
#'
#' # Construct a QFeatures object, with a 'peptideRaw' assay
#' pe_peptides <- make_pe(Silicosis_peptide, columns = ecols,   # columns is the abundance columns
#'                        expdesign = expDesign) # log2transform
#' pe_peptides # a QFeatures object, with a peptideRaw assay
make_pe <- function(Peptide, columns,
                    expdesign,
                    fnames,
                    assay_name = "peptideRaw",
                    log2transform = T){
  assertthat::assert_that(is.data.frame(Peptide), is.numeric(columns)|is.integer(columns)| is.character(columns),
                          is.data.frame(expdesign), is.character(assay_name), length(assay_name) == 1)
  if(is.numeric(columns)) columns = as.integer(columns)
  if(is.character(columns)){
    if(!all(columns %in% colnames(Peptide)))
      stop("columns should be the columns in ", deparse(substitute(Peptide)), "but ", columns[!columns %in% colnames(Peptide)],"do not exist.")
    columns = which(colnames(Peptide) %in% columns)
  }else if(is.integer(columns)){
    assert_that(all(columns %in% 1:nrow(Peptide)))
  }
  if (any(!apply(Peptide[, columns], 2, is.numeric))) {
    stop("Specified 'columns' should be numeric", "\nRun make_pe() with the appropriate columns as argument.",
         " Or run clean_character() to filter character values(like 'NA'), and transfer numeric in these columns",
         call. = FALSE)
  }

  if (any(!c("label", "condition", "replicate") %in%
          colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }

  if (tibble::is_tibble(Peptide))
    Peptide <- as.data.frame(Peptide)
  if (tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)

  rownames(expdesign) = expdesign$label

  raw = Peptide[,columns]
  expdesign <- mutate(expdesign, label = make.names(label), condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))) )
  if(any(is.na(matched))
     # && sum(matched,na.rm = T) != length(columns)
  ){
    matched <- match(make.names(expdesign$label),
                     make.names(colnames(raw)) )
  }
  if (any(is.na(matched))
      # && sum(matched,na.rm = T) != length(columns)
  ) {
    stop("Labels of the experimental design do not match ",
         "with column names in 'proteins_unique'", "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  colnames(Peptide)[columns] = colnames(raw)

  QF <- readQFeatures(table = Peptide, ecol = columns,
                      fnames = fnames,
                      name = assay_name)
  rowData(QF[[assay_name]])$nNonZero <- rowSums(assay(QF[[assay_name]]) > 0)
  # if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))

  QF[[assay_name]]@colData = as(expdesign, "DataFrame")
  colData(QF) = as(expdesign, "DataFrame")
  QF <- QFeatures::zeroIsNA(QF, assay_name)
  if(log2transform){ ## log2 Transform
    assay(QF[[1]]) = log2(assay(QF[[1]]))
  }
  return(QF)
}


#' Data.frame to QFeatures object conversion by parsing column names
#'
#' \code{make_pe_parse} creates a QFeatures object
#' based on a single data.frame. The experiment design is constructed
#' by colnames of expression assay.
#'
#' @inheritParams make_pe
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
#' @param remove_suffix Logical(1),
#' whether remove the suffix of expression columns.
#'
#' @return
#' An character(1) to name assay in the QFeatures object.
#' @export
#'
#' @examples
#' # Load a peptides
#' data(Silicosis_peptide)
#' ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = T)
#' # Construct a QFeatures object, with a 'peptideRaw' assay
#' pe_peptides <- make_pe_parse(Silicosis_peptide, columns = ecols,   # columns is the abundance columns
#'                              mode = "delim", sep = "_",
#'                              remove_prefix = T, log2transform = T)
#' pe_peptides
make_pe_parse <- function(Peptide,
                          columns,
                          fnames,
                          mode = c("char", "delim"),
                          chars = 1, sep = "_",
                          remove_prefix = T, remove_suffix = F,
                          assay_name = "peptideRaw",
                          log2transform = T){
  assertthat::assert_that(is.data.frame(Peptide), is.numeric(columns)|is.integer(columns)| is.character(columns),
                          is.character(assay_name), length(assay_name) == 1)
  mode <- match.arg(mode)

  if(is.numeric(columns)) columns = as.integer(columns)
  if(is.character(columns)){
    if(!all(columns %in% colnames(Peptide)))
      stop("columns should be the columns in ", deparse(substitute(Peptide)), "but ", columns[!columns %in% colnames(Peptide)],"do not exist.")
    columns = which(colnames(Peptide) %in% columns)
  }else if(is.integer(columns)){
    assert_that(all(columns %in% 1:nrow(Peptide)))
  }
  if (any(!apply(Peptide[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_pe_parse() with the appropriate columns as argument",
         call. = FALSE)
  }

  if (tibble::is_tibble(Peptide))
    Peptide <- as.data.frame(Peptide)

  if(remove_prefix){
    colnames(Peptide)[columns] <- delete_prefix(colnames(Peptide)[columns]) %>% make.names()
  }
  if(remove_suffix){
    colnames(raw) <- delete_suffix(colnames(raw)) %>% make.names()
  }
  if (mode == "char") {
    expdesign <- data.frame(label = colnames(Peptide)[columns], stringsAsFactors = FALSE) %>%
      mutate(condition = substr(label, 1, nchar(label) -
                                  chars), replicate = substr(label, nchar(label) +
                                                               1 - chars, nchar(label))) %>% unite(ID, condition,
                                                                                                   replicate, remove = FALSE)
  }
  if (mode == "delim") {
    # colnames(raw) = gsub(get_suffix(colnames(raw)),"", colnames(raw))
    expdesign <- data.frame(label = colnames(Peptide)[columns], stringsAsFactors = FALSE) %>%
      separate(label, c("condition", "replicate"), sep = sep,
               remove = FALSE, extra = "merge") %>% unite(ID,
                                                          condition, replicate, remove = FALSE)
  }

  rownames(expdesign) = expdesign$label
  Peptide[,columns] = apply(Peptide[,columns], 2, function(x){
    # temp1 <<- x
    # temp <<- as.integer(x)
    x[!(!grepl("[A-z]",x) & grepl("\\d",x))] = 0
    return(as.numeric(x))
  })

  QF <- readQFeatures(table = Peptide, ecol = columns,
                      fnames = fnames,
                      name = assay_name)
  rowData(QF[[assay_name]])$nNonZero <- rowSums(assay(QF[[assay_name]]) > 0)
  # if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))

  expdesign$replicate = as.integer(expdesign$replicate)
  QF[[assay_name]]@colData = as(expdesign, "DataFrame")
  colData(QF) = as(expdesign, "DataFrame")
  QF <- QFeatures::zeroIsNA(QF, assay_name)
  if(log2transform){ ## log2 Transform
    assay(QF[[1]]) = log2(assay(QF[[1]]))
  }
  return(QF)
}


#' Filter Qfeatures object on missing values or formula
#'
#' \code{filter_pe} filters a Qfeatures object based on missing values in assay
#' and provides formula.
#' The dataset is filtered for proteins that have a maximum of
#' 'thr' missing values in at least one condition; or a maximum of 'missnum' in all samples;
#' or a filter_formula
#'
#' @param pe Qfeatures,
#' Proteomics data (output from \code{\link{make_pe}()} or
#' \code{\link{make_pe_parse}()}).
#' @param assay_name Character(1),
#' the name of which assay to perform filter.
#' @inheritParams filter_se
#' @return A filtered Qfeatures object,
#' or a vector of keep rows if \code{return_keeprows} is TRUE
#'
#' @examples
#' # Load example peptide data
#' data(Silicosis_peptide)
#' ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = TRUE)
#'
#' # Construct a QFeatures object, with a 'peptideRaw' assay
#' pe_peptides <- make_pe_parse(Silicosis_peptide, columns = ecols, remove_prefix = TRUE, log2transform = TRUE,mode = "delim")
#' filt_pe <- filter_pe(pe_peptides, thr = 1,fraction = 0.4, filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" )
#'
#' @export
#' @importFrom rlang is_formula
filter_pe <- function(pe,
                      thr = NULL,
                      missnum = NULL,
                      fraction = NULL,
                      filter_formula = NULL,
                      assay_name = "peptideRaw",
                      return_keeprows = F){
  if(is.character(filter_formula) & length(filter_formula) == 1) filter_formula = as.formula(filter_formula)
  assertthat::assert_that(class(pe) == "QFeatures", is.null(thr)|is.numeric(thr), is.null(missnum)|is.numeric(missnum),
                          is.character(assay_name), length(assay_name) == 1, rlang::is_formula(filter_formula),
                          is.null(fraction)|is.numeric(fraction), is.null(fraction)|length(fraction) == 1)

  if(!assay_name %in% names(pe))
    stop("'assay_name' should be one of exist assay in ", deparse(substitute(pe)),": ", paste(names(pe),collapse = ", "))
    # stop("'assay_name' do not exist in pe, please validate")

  se <- pe[["peptideRaw"]]
  keep_row = filter_se(se, thr = thr, missnum = missnum, fraction = fraction,
                       filter_formula = filter_formula,return_keeprows = T)
  if(return_keeprows) return(keep_row)
  pe = pe[keep_row,]
  return(pe)
}



# aggregation functions
## modified from QFeatures:::reduceDataFrame, add reserve to trim the row DF.
#' @import  QFeatures
#' @importFrom S4Vectors split DataFrame isEmpty
reducedataframe <- function (x, k, count = FALSE, simplify = TRUE, drop = FALSE, reserve)
{
  res <- S4Vectors::split(x, k)
  lens <- unname(lengths(res))
  colnames(x)

  if (simplify | drop)
    invars <- invariant_cols2(res)
  res <- S4Vectors::DataFrame(res)

  if(!is.null(reserve)){
    invars <- c(invars,which(colnames(res) %in% reserve)) %>% unique
  }
  if (simplify) {
    for (i in invars) res[[i]] <- unname(sapply(res[[i]],
                                                "[[", 1))
  }

  if (drop)
    res <- res[, invars, drop = FALSE]

  if (count)
    res[[".n"]] <- lens

  res
}

## from QFeatures:::invariant_cols2
invariant_cols2 <- function (x)
{
  res <- rep(NA, length(x[[1]]))
  for (i in seq_along(res)) {
    res[i] <- all(lengths(lapply(x[, i], unique)) == 1)
  }
  which(res)
}

## modified from QFeatures:::.aggregateQFeatures, add reserve argument for reducedataframe
#' @importFrom MsCoreUtils aggregate_by_vector robustSummary colCounts
.aggregateSE <- function(object, fcol, fun, reserve,...) {
  if (missing(fcol))
    stop("'fcol' is required.")
  m <- assay(object, 1)
  rd <- rowData(object)
  if (!fcol %in% names(rd))
    stop("'fcol' not found in the assay's rowData.")
  groupBy <- rd[[fcol]]

  ## Store class of assay i in case it is not a Summarized experiment so that
  ## the aggregated assay can be reverted to that class
  .class <- class(object)

  ## Message about NA values is quant/row data
  has_na <- character()
  if (anyNA(m))
    has_na <- c(has_na, "quantitative")
  if (anyNA(rd, recursive = TRUE))
    has_na <- c(has_na, "row")
  if (length(has_na)) {
    msg <- paste(paste("Your", paste(has_na, collapse = " and "),
                       " data contain missing values."),
                 "Please read the relevant section(s) in the",
                 "aggregateFeatures manual page regarding the",
                 "effects of missing values on data aggregation.")
    message(paste(strwrap(msg), collapse = "\n"))
  }

  aggregated_assay <- MsCoreUtils::aggregate_by_vector(m, groupBy, fun, ...)
  aggcount_assay <- MsCoreUtils::aggregate_by_vector(m, groupBy, colCounts)
  print("begin reducedataframe")
  aggregated_rowdata <- reducedataframe(rd, rd[[fcol]],
                                        simplify = T,
                                        drop = T,
                                        count = TRUE,
                                        reserve = reserve)

  print("reducedataframe finiched")
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(assay = aggregated_assay,
                                                            aggcounts = aggcount_assay),
                             rowData = aggregated_rowdata[rownames(aggregated_assay), ])
  ## If the input objects weren't SummarizedExperiments, then try to
  ## convert the merged assay into that class. If the conversion
  ## fails, keep the SummarizedExperiment, otherwise use the
  ## converted object .
  if (.class != "SummarizedExperiment")
    se <- tryCatch(as(se, .class),
                   error = function(e) se)

  return(se)
}


## distribute razor peptides to proteingroups
Peptide_distribution <- function(pe_norm, i = "peptideNorm", fcol = "Proteins"){
  rd <- rowData(pe_norm[[i]]) %>% as.data.frame()
  allpgs <- rd[,fcol]

  ## get the smallestUnique PGs
  SU_PGs <- smallestUniqueGroups(allpgs)
  SU_PGs_uniquecounts <- allpgs[allpgs %in% SU_PGs] %>% table ## The unique peptide number in each SU PGs
  SU_PGs_uniquecounts %>% class
  SU_PGs_list1 <- strsplit(SU_PGs,";",fixed = F)
  names(SU_PGs_list1) <- SU_PGs

  SU_PGs_df <- lapply(1:length(SU_PGs_list1), function(x){
    return(data.frame(suPG = names(SU_PGs_list1)[x], proteinsID = SU_PGs_list1[[x]]))
  }) %>% do.call(rbind,.)
  SU_PGs_df$uniquecounts = SU_PGs_uniquecounts[match(SU_PGs_df$suPG,names(SU_PGs_uniquecounts))]
  SU_PGs_df = arrange(SU_PGs_df,desc(uniquecounts))

  rd$peptide.type <- ifelse(rd[,fcol] %in% SU_PGs, "unique", "razor")
  nonSU_PGs <- rd %>% filter(peptide.type == "razor") %>% .[,fcol] %>% unique() ## non-smallestUnique PGs

  ({
    nonSU_PGs_list3 = nonSU_PGs %>% sapply(., function(x){
      x = strsplit(x,";")[[1]]
      temp = subset(SU_PGs_df,proteinsID %in% x)
      temp = temp[temp$uniquecounts == temp$uniquecounts[1],]
      razorPG = temp$suPG[match(x,temp$proteinsID) %>% na.omit %>% .[1]]
      return(razorPG)
    })
  }) %>% system.time # to calulate spend time
  # class(nonSU_PGs_list3)

  rd$smallestProteingroups <- ifelse(rd$peptide.type=="unique" , rd[,fcol], nonSU_PGs_list3[match(rd[,fcol],names(nonSU_PGs_list3))])
  rowData(pe_norm[[i]]) <- rd
  return(pe_norm)
}


## The function from msqrob2, design the smallestUniqueGroups. See example below.
#> DEP2:::smallestUniqueGroups(c("A;B;C","D","B;C","A;B","B;C;D"))
#> [1] "D"   "B;C" "A;B"
smallestUniqueGroups <- function(proteins,
                                 split = ";") {
  b <- strsplit(x = as.character(unique(proteins)), split = split, fixed = TRUE)
  ## b is a list store proteins names(or ids)
  included <- vector()

  j <- 1 # j control how many proteins in a single names vector(protein group)
  while (length(b) != 0) {
    # To exctract the names vectors contain j proteins from list b.
    included <- c(
      included,
      vapply(
        b[vapply(b, length, integer(1)) == j],
        function(x) paste(x, collapse = split),
        character(1)
      )
    )
    a <- unlist(b[vapply(b, length, integer(1)) == j])
    # To exctract the names vectors contain more than j proteins from list b.
    b <- b[vapply(b, length, integer(1)) > j] # remove names with <= j proteins from list b

    if (length(b) != 0) {
      sel <- vector()
      for (i in seq_len(length(b))) {
        sel[i] <- !any(b[[i]] %in% a) # if names contains a protein overlapped in included F, else T.
      }
      b <- b[sel] # remove long pgs with any protein overlapped with short PGs. Make long pgs overlap with short PGs.
      j <- j + 1 # increase j
    }
  }

  included <- unlist(included)
  return(included)
}


aggregateFeatures = function(object, i, fcol, name = "newAssay",
                           fun = MsCoreUtils::robustSummary, ...) {
  if (S4Vectors::isEmpty(object))
    return(object)
  if (name %in% names(object))
    stop("There's already an assay named '", name, "'.")
  if (missing(i))
    i <- QFeatures:::main_assay(object)
  print("begin assay aggregate")
  ## Create the aggregated assay
  aggAssay <- .aggregateSE(object[[i]], fcol, fun, ...)
  print("assay aggregate finished")
  ## Add the assay to the QFeatures object
  object <- QFeatures::addAssay(object,
                                aggAssay,
                                name = name)
  ## Link the input assay to the aggregated assay
  QFeatures::addAssayLink(object,
               from = i,
               to  = name,
               varFrom = fcol,
               varTo = fcol)
}


#' Normalize a QFeatures object
#'
#' Normalize a QFeatures object though [QFeatures::normalize] function
#'
#' @param pe A QFeature object
#' @param method 	Character(1), normalisation method, one of "diff.median", "quantiles", "quantiles.robust" or "vsn".
#' See \link[QFeatures]{normalize}.
#' @param i 	A numeric vector or a character vector giving the index or the name, respectively, of the assay(s) to be processed.
#' @param name  Character(1) naming the new normalized assay name.
#'
#' @importFrom QFeatures normalize
#'
#' @examples
#' \dontrun{
#' # Load example peptide data
#' data(Silicosis_peptide)
#' ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = T)
#'
#' # Construct QFeature object
#' pe_peptides <- make_pe_parse(Silicosis_peptide, columns = ecols, remove_prefix = T, log2transform = T,mode = "delim")
#'
#' # Filter, imputation
#' filt_pe <- filter_pe(pe_peptides, thr = 1,fraction = 0.4, filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" )
#' imp_pe <- QFeatures::addAssay(filt_pe, DEP2::impute(filt_pe[["peptideRaw"]], fun = "MinDet"), name = "peptideImp")
#'
#' # Normalization
#' norm_pe <- normalize_pe(imp_pe,method = "quantiles", i = "peptideImp", name = "peptideNorm")
#' }
#' @export
normalize_pe <- function(pe, method = c("diff.median", "quantiles", "quantiles.robust" ,"vsn"), i = "peptideRaw", name = "peptideNorm"){
  method = match.arg(method)
  if(is.character(i) && !i %in% names(pe)) stop("'i' should be one of exist assay in ", deparse(substitute(pe)),": ", paste(names(pe),collapse = ", "))
  if(is.numeric(i) | is.integer(i)) assertthat::assert_that(i %in% 1:length(pe))
  if(method == "vsn"){
    assay(pe[[i]]) = 2^assay(pe[[i]])
    pe_norm <- QFeatures::normalize(pe,
                                    i = i,
                                    name = name,
                                    method = method) ## should be MsCoreUtils::normalizeMethods
    assay(pe_norm[[i]]) = log2(assay(pe_norm[[i]]))
  }else if(method %in% c("quantiles", "diff.median")){
    pe_norm <- QFeatures::normalize(pe,
                                    i = i,
                                    name = name,
                                    method = method)
  }else if(method %in% c("quantiles.robust")){
    if(anyNA(assay(pe[[i]]))) stop("quantiles.robust normalize can work for data with missing values (NA).")
    pe_norm <- QFeatures::normalize(pe,
                                    i = i,
                                    name = name,
                                    method = method)
  }
  return(pe_norm)
}


#' Summarize peptide quantity to protein quantity
#'
#' Aggregation peptide to protein vis 'aggrefun'.
#'
#' @param pe A QFeatures object,
#' contains the normalized peptide assay
#' @param aggrefun 	A function used for quantitative feature aggregation.
#' It can be a character in "RobustSummary","medianPolish","totalMean" or other function.
#' Details see \code{\link[QFeatures]{aggregateFeatures}}
#' @param aggregate_Peptide_Type Character in "Unique + Razor" or "Unique".
#' Use what kind of peptides to summarise proteins. If choose "Unique", return output just save unique peptides in smallest proteingroups.
#' @param fcol Character(1), defining how to summarise the features. Exist in \code{rowData(pe)}.
#' @param peptide_assay_name Character(1), the name of aggregation result assay
#' @param reserve Character, the column(s) which will reserve after aggregate, such as the columns store protein information can.
#'
#' @return
#' A QFeatures object with a new protein aggregation assay.
#' @examples
#' \dontrun{
#' # Load example peptide data
#' data(Silicosis_peptide)
#' ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = T)
#'
#' # construct QFeatures object
#' pe_peptides <- make_pe_parse(Silicosis_peptide, columns = ecols, remove_prefix = T, log2transform = T,mode = "delim")
#' filt_pe <- filter_pe(pe_peptides, thr = 1,fraction = 0.4, filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" )
#' imp_pe <- QFeatures::addAssay(filt_pe, DEP2::impute(filt_pe[["peptideRaw"]], fun = "MinDet"), name = "peptideImp")
#' norm_pe <- DEP2:::normalize_pe(imp_pe,method = "quantiles", i = "peptideImp", name = "peptideNorm")
#'
#' # Summarize peptide value to protein quantity
#' protein_pe <- DEP2::aggregate_pe(norm_pe, fcol = "Proteins", peptide_assay_name = "peptideNorm")
#' class(protein_pe)
#'
#' }
#' @export
#'
aggregate_pe <- function(pe, aggrefun = c("RobustSummary","medianPolish","totalMean"), aggregate_Peptide_Type = c("Unique + Razor", "Unique"),
                         fcol, peptide_assay_name = "peptideNorm", reserve = "Gene.names"
){
  aggregate_Peptide_Type <- match.arg(aggregate_Peptide_Type)
  assertthat::assert_that(class(pe) == "QFeatures", is.function(aggrefun)|is.character(aggrefun),
                          is.character(aggregate_Peptide_Type), length(aggregate_Peptide_Type) == 1, is.character(aggregate_Peptide_Type))
  if(class(aggrefun) == "character"){
    aggrefun <- match.arg(aggrefun,choices = c("RobustSummary","medianPolish","totalMean"))
    aggrefun = switch(aggrefun,
                      totalSum = base::colSums,
                      totalMean = base::colMeans,
                      medianPolish = MsCoreUtils::medianPolish,
                      RobustSummary = MsCoreUtils::robustSummary)
  }

  if(!peptide_assay_name %in% names(pe))
    stop("'peptide_assay_name' should be one of exist assay in ", deparse(substitute(pe)),": ", paste(names(pe),collapse = ", "))

  if(aggregate_Peptide_Type == "Unique"){
    ## just use the features in smallestUniqueGroups.
    fil_formula <- as.formula( paste0("~",fcol," %in% smallestUniqueGroups(rowData(pe[['",peptide_assay_name ,"']])$",fcol,")") )
    pe <- filter_pe(pe, filter_formula = fil_formula, assay_name =  peptide_assay_name )
    print("aggregate by uniques peptides, filterFeatures finished")
    rowData(pe)$smallestProteingroups <- rowData(pe)[,fcol]
    protein <- suppressWarnings(
      {aggregateFeatures(object = pe ,
                         i = peptide_assay_name, fcol = fcol,
                         name = "protein",
                         fun = aggrefun,
                         na.rm = T,
                         reserve = reserve)}
      # QFeatures::aggregateFeatures(object = pe ,
      #                              i = peptide_assay_name, fcol = fcol,
      #                              name = "protein",
      #                              fun = aggrefun,
      #                              na.rm = T)
    )


  }else if(aggregate_Peptide_Type == "Unique + Razor"){
    print("aggregate by Unique + Razor peptides")
    pe <- Peptide_distribution(pe, i = peptide_assay_name,fcol = fcol)
    print("peptides distribution finished")
    protein <- suppressWarnings({
      aggregateFeatures(object = pe ,
                        i = peptide_assay_name, fcol = "smallestProteingroups",
                        name = "protein",
                        fun = aggrefun,
                        na.rm = T,
                        reserve = reserve)
      # QFeatures::aggregateFeatures(object = pe ,
      #                              i = peptide_assay_name, fcol = "smallestProteingroups",
      #                              name = "protein",
      #                              fun = aggrefun,
      #                              na.rm = T)
    } ) ## function from QFeatures
  }
  print("aggregation finished")
  colData(protein[["protein"]]) = colData(protein)
  return(protein)
}



#' Extract the proteins SummarizedExperiment object from a QFeatures container
#'
#' This function accept a result from \link{aggregate_pe}(),
#' tidy features identifiers and return a SE object for following analyze
#'
#' @param pe_aggregated A QFeatures object output from \code{aggregate_pe()}, which contain a "protein" quantative assay.
#' @param names The column of gene names, which serve as feature identifier and is transmitted to \link{make_unique}
#' @param ids   The column of protein ID, transmitted to \link{make_unique}. aggregate_pe automatically generate a "smallestProteingroups" column
#' to store ids in proteingroups.
#' @param delim Character. The separator in names of ids.
#'
#' @return
#' A SummarizedExperiment object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example peptide data
#' data(Silicosis_peptide)
#' ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = T)
#'
#' # construct QFeatures object
#' pe_peptides <- make_pe_parse(Silicosis_peptide, columns = ecols, remove_prefix = T, log2transform = T,mode = "delim")
#' filt_pe <- filter_pe(pe_peptides, thr = 1,fraction = 0.4, filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" )
#' imp_pe <- QFeatures::addAssay(filt_pe, DEP2::impute(filt_pe[["peptideRaw"]], fun = "MinDet"), name = "peptideImp")
#' norm_pe <- DEP2:::normalize_pe(imp_pe,method = "quantiles", i = "peptideImp", name = "peptideNorm")
#'
#' # Summarize peptide value to protein quantity
#' protein_pe <- DEP2::aggregate_pe(norm_pe, fcol = "Proteins", peptide_assay_name = "peptideNorm")
#' class(protein_pe)
#'
#' # Construct a SE object
#' se <- pe2se(protein_pe)
#' class(se)
#' }
pe2se <- function(pe_aggregated, names = "Gene.names", ids = "smallestProteingroups", delim = ";"){
  SE_pep <- pe_aggregated[["protein"]]

  ## make names for features
  rd <- DEP2::make_unique(rowData(SE_pep) %>% as.data.frame(),
                          names = names,           # The gene names
                          ids = ids,  # The protien ID translated to Proteingroups.
                          delim = ";")
  rownames(rd) <- rd$name
  rowData(SE_pep) <- rd
  rownames(SE_pep) <- rd$name

  colData(SE_pep)$replicate <- as.integer(colData(SE_pep)$replicate)
  return(SE_pep)
}

