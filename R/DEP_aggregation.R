

#' Data.frame to QFeatures object conversion using an experimental design
#'
#' \code{make_pe} creates a QFeatures object
#' based on two data.frames: the peptide table and experimental design.
#'
#' @param Peptide
#' Data.frame,
#' Peptide table with its belonging proteingroup.
#' @param ecols A numeric indicating the indices of the columns
#' to be used as expression values. Can also be a character indicating the names of the columns.
#' @param protein_col An optional character(1) or numeric(1) indicating the column of
#' peptides' proteins/proteingroups information.
#' @param expdesign A data.frame.
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' @param assay_name
#' An character(1) to name assay in the QFeatures object.
#' @param log2transform Logical(1), whether log2 transform the assay, default TRUE.
#' @return
#' A QFeatures object with a peptide assay.
#' @export
#'
#' @examples
make_pe <- function(Peptide, ecols, protein_col,
                    expdesign,
                    assay_name = "peptideRaw",
                    log2transform = T){
  if(is.numeric(ecols)) ecols = as.integer(ecols)
  assertthat::assert_that(is.data.frame(Peptide), is.integer(ecols),
                          is.data.frame(expdesign), is.character(assay_name), length(assay_name) == 1)

  if (any(!apply(Peptide[, ecols], 2, is.numeric))) {
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

  raw = Peptide[,ecols]
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
  colnames(Peptide)[ecols] = colnames(raw)

  QF <- readQFeatures(table = Peptide, ecol = ecols,
                      fnames = protein_col,
                      name = assay_name)
  rowData(QF[[assay_name]])$nNonZero <- rowSums(assay(QF[[assay_name]]) > 0)
  # if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))

  QF[[assay_name]]@colData = as(expdesign, "DataFrame")
  colData(QF) = as(expdesign, "DataFrame")
  QF <- zeroIsNA(QF, assay_name)
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
#'
#' @return
#' An character(1) to name assay in the QFeatures object.
#' @export
#'
#' @examples
make_pe_parse <- function(Peptide,
                          ecols,
                          protein_col,
                          mode = c("char", "delim"),
                          chars = 1, sep = "_", remove_prefix = T,
                          assay_name = "peptideRaw",
                          log2transform = T){
  mode <- match.arg(mode)
  if (any(!apply(Peptide[, ecols], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_pe_parse() with the appropriate columns as argument",
         call. = FALSE)
  }

  if (tibble::is_tibble(Peptide))
    Peptide <- as.data.frame(Peptide)

  if(remove_prefix){
    colnames(Peptide)[ecols] <- DEP2:::delete_prefix(colnames(Peptide)[ecols]) %>% make.names()
  }
  if (mode == "char") {
    expdesign <- data.frame(label = colnames(Peptide)[ecols], stringsAsFactors = FALSE) %>%
      mutate(condition = substr(label, 1, nchar(label) -
                                  chars), replicate = substr(label, nchar(label) +
                                                               1 - chars, nchar(label))) %>% unite(ID, condition,
                                                                                                   replicate, remove = FALSE)
  }
  if (mode == "delim") {
    # colnames(raw) = gsub(get_suffix(colnames(raw)),"", colnames(raw))
    expdesign <- data.frame(label = colnames(Peptide)[ecols], stringsAsFactors = FALSE) %>%
      separate(label, c("condition", "replicate"), sep = sep,
               remove = FALSE, extra = "merge") %>% unite(ID,
                                                          condition, replicate, remove = FALSE)
  }

  rownames(expdesign) = expdesign$label
  Peptide[,ecols] = apply(Peptide[,ecols], 2, function(x){
    # temp1 <<- x
    # temp <<- as.integer(x)
    x[!(!grepl("[A-z]",x) & grepl("\\d",x))] = 0
    return(as.numeric(x))
  })

  QF <- readQFeatures(table = Peptide, ecol = ecols,
                      fnames = protein_col,
                      name = assay_name)
  rowData(QF[[assay_name]])$nNonZero <- rowSums(assay(QF[[assay_name]]) > 0)
  # if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))

  QF[[assay_name]]@colData = as(expdesign, "DataFrame")
  colData(QF) = as(expdesign, "DataFrame")
  QF <- zeroIsNA(QF, assay_name)
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
#' @examples
#'
#' @export
#' @importFrom rlang is_formula
filter_pe <- function(pe,
                      thr = NULL,
                      missnum = NULL,
                      filter_formula = NULL,
                      assay_name = "peptideRaw",
                      return_keeprows = F){
  if(is.character(filter_formula) & length(filter_formula) == 1) filter_formula = as.formula(filter_formula)
  assertthat::assert_that(class(pe) == "QFeatures", is.null(thr)|is.numeric(thr), is.null(missnum)|is.numeric(missnum),
                          is.character(assay_name), length(assay_name) == 1, rlang::is_formula(filter_formula))

  if(!assay_name %in% names(pe))
    stop("'assay_name' should be one of exist assay in ", deparse(substitute(pe)),": ", paste(names(pe),collapse = ", "))
    # stop("'assay_name' do not exist in pe, please validate")

  se <- pe[["peptideRaw"]]
  keep_row = filter_se(se, thr = thr, missnum = missnum, filter_formula = filter_formula,return_keeprows = T)
  if(return_keeprows) return(keep_row)
  pe = pe[keep_row,]
  return(pe)
}



# aggregation functions
## modified from QFeatures:::reduceDataFrame, add reserve.
#' @import  QFeatures
#' @importFrom S4Vectors split
reducedataframe <- function (x, k, count = FALSE, simplify = TRUE, drop = FALSE, reserve)
{
  # x_save <<- x
  # k_save <<- k
  # count_save <<- count
  # simplify_save <<- simplify
  # drop_save <<- drop
  # reserve_save <<- reserve
  # cat("aaaa1")
  res <- S4Vectors::split(x, k)
  # res_save <<- res
  lens <- unname(lengths(res))
  # lens_save1 <<- lens
  colnames(x)
  # cat("aa")
  if (simplify | drop)
    invars <- QFeatures:::invariant_cols2(res)
  res <- DataFrame(res)
  # res_save1 <<- res
  # cat("bb")
  if(!is.null(reserve)){
    invars <- c(invars,which(colnames(res) %in% reserve)) %>% unique
  }
  if (simplify) {
    for (i in invars) res[[i]] <- unname(sapply(res[[i]],
                                                "[[", 1))
  }
  # res_save2 <<- res
  # cat("cc")
  if (drop)
    res <- res[, invars, drop = FALSE]
  # cat("dd")
  # res_save3 <<- res
  # lens_save2 <<- lens
  if (count)
    res[[".n"]] <- lens
  # cat("ee")
  res
}



## modified from QFeatures:::.aggregateQFeatures, add reserve for reducedataframe
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
  aggregated_rowdata <- DEP2:::reducedataframe(rd, rd[[fcol]],
                                        simplify = T,
                                        drop = T,
                                        count = TRUE,
                                        reserve = reserve)

  print("reducedataframe finiched")
  se <- SummarizedExperiment(assays = SimpleList(assay = aggregated_assay,
                                                 aggcounts = aggcount_assay),
                             rowData = aggregated_rowdata[rownames(aggregated_assay), ])
  ## If the input objects weren't SummarizedExperiments, then try to
  ## convert the merged assay into that class. If the conversion
  ## fails, keep the SummarizedExperiment, otherwise use the
  ## converted object (see issue #78).
  if (.class != "SummarizedExperiment")
    se <- tryCatch(as(se, .class),
                   error = function(e) se)

  return(se)
}
# environment(.aggregateSE) = asNamespace("QFeatures")


## distribute razor peptides
Peptide_distribution <- function(pe_norm, i = "peptideNorm", fcol = "Proteins"){
  rd <- rowData(pe_norm[[i]]) %>% as.data.frame()
  allpgs <- rd[,fcol]

  ## get the smallestUnique PGs
  SU_PGs <- smallestUniqueGroups(allpgs)
  SU_PGs_uniquecounts <- allpgs[allpgs %in% SU_PGs] %>% table ## 计算每个最终 PG 的unique数量
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
    })}) %>% system.time
  class(nonSU_PGs_list3)

  rd$smallestProteingroups <- ifelse(rd$peptide.type=="unique" , rd[,fcol], nonSU_PGs_list3[match(rd[,fcol],names(nonSU_PGs_list3))])
  rowData(pe_norm[[i]]) <- rd
  return(pe_norm)
}


## The function from msqrob2, design the smallestUniqueGroups
smallestUniqueGroups <- function(proteins,
                                 split = ";") {
  b <- strsplit(x = as.character(unique(proteins)), split = split, fixed = TRUE)

  included <- vector()

  j <- 1
  while (length(b) != 0) {
    #提出b这个list中length是j的蛋白名,included
    included <- c(
      included,
      vapply(
        b[vapply(b, length, integer(1)) == j],
        function(x) paste(x, collapse = split),
        character(1)
      )
    )
    a <- unlist(b[vapply(b, length, integer(1)) == j])
    ##提出b这个list中length > j的蛋白名
    b <- b[vapply(b, length, integer(1)) > j]

    if (length(b) != 0) {
      sel <- vector()
      for (i in seq_len(length(b))) {
        sel[i] <- !any(b[[i]] %in% a) #length是1的蛋白名,即included不在list中length > 1的蛋白名中为TRUE，在的话为FALSE
      }
      b <- b[sel] #把length是1的蛋白名,即included不在list中length > 1的蛋白名中为TRUE，提出来，重新赋予b
      j <- j + 1 #j=2再返回while循环
    }
  }

  included <- unlist(included)
  return(included)
}


aggregateFeatures = function(object, i, fcol, name = "newAssay",
                           fun = MsCoreUtils::robustSummary, ...) {
  if (isEmpty(object))
    return(object)
  if (name %in% names(object))
    stop("There's already an assay named '", name, "'.")
  if (missing(i))
    i <- QFeatures:::main_assay(object)
  print("begin assay aggregate")
  ## Create the aggregated assay
  aggAssay <- DEP2:::.aggregateSE(object[[i]], fcol, fun, ...)
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
#' @param i 	A numeric vector or a character vector giving the index or the name, respectively, of the assay(s) to be processed.
#' @param name  Character(1) naming the new normalized assay name.
#' @importFrom QFeatures normalize
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

# aggregate_pe(pe_save, aggrefun = "RobustSummary", aggregate_Peptide_Type = "Unique + Razor",fcol = "Proteins",
#              peptide_assay_name = "peptideNorm", reserve = "Gene.names")


#' Summarize peptide quantity to protein quantity
#'
#' Aggregation peptide to protein vis 'aggrefun'.
#'
#' @param pe A QFeatures object,
#' contains the normalized peptide assay
#' @param aggrefun 	A function used for quantitative feature aggregation.
#' It can be a character in "RobustSummary","medianPolish","totalMean" or a function.
#' Details see \code{\link{QFeatures::aggregateFeatures}}
#' @param aggregate_Peptide_Type Character in "Unique + Razor" or "unique".
#' Use what kind of peptides to summarise proteins.
#' @param fcol Character(1), defining how to summarise the features. Exist in \code{rowData(pe)}.
#' @param peptide_assay_name Character(1), the name of aggregation result assay
#' @param reserve Character, the column(s) which will reserve after aggregate.
#'
#' @return
#' a QFeatures object with a new protein aggregation assay.
#' @export
#'
#' @examples
aggregate_pe <- function(pe, aggrefun = c("RobustSummary","medianPolish","totalMean"), aggregate_Peptide_Type = c("Unique + Razor", "unique"),
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

  if(aggregate_Peptide_Type == "unique"){
    fil_formula <- as.formula( paste0("~",fcol," %in% smallestUniqueGroups(rowData(pe[['",peptide_assay_name ,"']])$",fcol,")") )
    pe <- filter_pe(pe, fil_formula, assay_name =  peptide_assay_name )
    print("aggregate by uniques peptides, filterFeatures finished")

    protein = suppressWarnings({DEP2:::aggregateFeatures(object = pe ,
                                                  i = peptide_assay_name, fcol = "smallestProteingroups",
                                                  name = "protein",
                                                  fun = aggrefun,
                                                  na.rm = T,
                                                  reserve = reserve)} )


  }else if(aggregate_Peptide_Type == "Unique + Razor"){
    print("aggregate by Unique + Razor peptides")
    pe <- Peptide_distribution(pe, i = peptide_assay_name,fcol = fcol)
    print("peptides distribution finished")
    protein = suppressWarnings({DEP2:::aggregateFeatures(object = pe ,
                                                         i = peptide_assay_name, fcol = "smallestProteingroups",
                                                         name = "protein",
                                                         fun = aggrefun,
                                                         na.rm = T,
                                                         reserve = reserve)} ) ## function from QFeatures
  }
  print("aggregation finished")
  colData(protein[["protein"]]) = colData(protein)
  return(protein)
}



