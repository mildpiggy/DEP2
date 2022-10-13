

# extract significant candidate in giving contrast. If contrast NULL, extract all significant candidate. Return a SE.
#' Extract significant candidates from SummarizedExperiment or DEGdata
#'
#' @param object a SummarizedExperiment or DEGdata object.
#' @param contrasts NULL or contrasts in object.
#' @param thresholdmethod NULL or 'intersect' or 'curve'. The thresholdmethod to decide significant.
#' If is NULL, used existing rejections. Otherwise filter new rejections via
#' \code{\link{add_rejections}}().
#' @param curvature Numeric(1), Sets the curvature for the curve cutoff lines
#' @param x0_fold Numeric(1), decide the x0 ratio to the standard deviations of L2FC. The x0 usually is set to 1(medium confidence) or 2(high confidence) standard deviations.
#' @param return_type One of "subset", "table", "names". "subset" return a subset object, "table" return a result data.frame,
#' names return name vector of object.
#'
#' @return
#' A object in the same class of input or a data.frame or a names vector, designed by return_type.
#'
#' @export
#'
#' @examples
get_signicant <- function(object,
                          contrasts = NULL,
                          thresholdmethod = NULL,
                          curvature = 1,
                          x0_fold = 2,
                          return_type = c("subset", "table", "names")
                          ){
  row_data = rowData(object)
  return_type = match.arg(return_type)
  # trend = match.arg(trend)
  if(is.null(thresholdmethod)){
    ## check significant cols
    if(length(grep("_significant$",colnames(row_data))) < 1){
      stop(paste0("'significant' column is necessary when thresholdmethod = NULL, and is not presented in '",
                  deparse(substitute(object)), "'.\nRun add_rejections() to obtain the required column,
                  or set a significant thresholdmethod to obtain rejections."))
    }
    de = object
  }else{
    ## check test results
    if(length(grep("_diff$",colnames(row_data))) < 1){
      if(class(object) == "SummarizedExperiment") frontfun = "test_diff"
      stop(paste0("Differentail test result is necessary for classify significant candidate via 'thresholdmethod',
                  and is not presented in '",
                  deparse(substitute(object)), "'.\nRun a ", frontfun, "() to obtain the required column,
                  or set a significant thresholdmethod."))
    }
    de = add_rejections(object,
                        # contrasts = contrasts,
                        thresholdmethod = thresholdmethod,
                        curvature = 1,
                        x0_fold = 2)
    de_save <<- de
  }

  row_data = rowData(de)
  if(is.null(contrasts)){
    r_ind = row_data[,"significant"]
  }else{
    exist_contrasts <- get_contrast(de)
    ## check input contrasts
    if(any(!contrasts %in% exist_contrasts)){
      if(all(!contrasts %in% exist_contrasts)){
        stop("Input contrasts: ", paste0(contrasts, collapse = ", ")," don't exists. Contrasts should be: ",
             paste0(exist_contrasts, collapse = ", "))
      }else{
        warning("Input contrasts: ", paste0(contrasts[which(!contrasts %in% exist_contrasts)], collapse = ", ")," don't exists. Only perform on: ",
                paste0(intersect(contrasts, exist_contrasts), collapse = ", "))
        contrasts = intersect(contrasts, exist_contrasts)
      }
    }
    sig_ind <- row_data[ , paste(contrasts, "_significant", sep = ""),drop = F]
    sig_ind_save <<- sig_ind
    if(class(sig_ind) == "data.frame" || class(sig_ind) == "DFrame"){
      r_ind = apply(sig_ind, 1, any)

    }
  }

  filtered <- de[r_ind , ]
  row_data_filtered = rowData(filtered)
  return_val = switch(return_type,
                      "subset" = filtered,
                      "table" = get_results(filtered),
                      "names" = rownames(row_data_filtered)
                      )
}



# get_signicant(dds_res, contrasts = get_contrast(dds_res)[1], thresholdmethod = "curve") -> temp
# temp


# get_contrast from dep
#' Get contrast(s) from SummarizedExperiment or DEGdata
#'
#' @param object A SummarizedExperiment or DEGdata object output from \link{test_diff}()
#'
#' @return
#' Character, the existing contrast(s) in object, which is introduced in the test_diff step.
#'
#' @export
#'
#' @examples
get_contrast = function(object){
  row_data = rowData(object)
  if (length(grep("_diff$", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '",
                deparse(substitute(object)), "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  colnames(row_data) %>% grep("_diff$", .,value = TRUE) %>% gsub("_diff$", "",.)
}


# from DEP::exclude_deps
# Internal function to exclude differentially expressed candidates from
# particular contrasts. Accept both SE or DEGdata.
exclude_contrasts <- function(object, contrasts) {
  assertthat::assert_that(inherits(object, "SummarizedExperiment"))

  row_data <- rowData(object, use.names = FALSE)
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         objectarse(substitute(object)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }

  if(is.null(contrasts)) {
    filtered <- object
  } else {
    contrasts_colnames <- paste0(contrasts, "_significant")
    matches <- match(contrasts_colnames, colnames(row_data))

    if (any(is.na(matches))) {
      valid_cntrsts <- row_data %>%
        data.frame() %>%
        select(ends_with("_diff")) %>%
        colnames(.) %>%
        gsub("_diff", "", .)
      valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                  paste0(valid_cntrsts, collapse = "', '"),
                                  "'")
      stop("The contrast(s) is/are not valid, ",
           "please run `exclude_contrasts()` with a valid contrast as argument\n",
           valid_cntrsts_msg,
           call. = FALSE)
    }

    if(length(matches) == 1) {
      filtered <- object[!row_data[,matches],]
    } else {
      filtered <- object[!apply(row_data[,matches], 1, any)]
    }
  }
  return(filtered)
}

# from DEP::select_deps
# Internal function to select differentially expressed candidates from
# particular contrasts. Accept both SE or DEGdata.
select_contrasts <- function(object, contrasts) {
  assertthat::assert_that(inherits(object, "SummarizedExperiment"))

  row_data <- rowData(object, use.names = FALSE)
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         objectarse(substitute(object)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }

  if(is.null(contrasts)) {
    filtered <- object
  } else {
    contrasts_colnames <- paste0(contrasts, "_significant")
    matches <- match(contrasts_colnames, colnames(row_data))

    if (any(is.na(matches))) {
      valid_cntrsts <- row_data %>%
        data.frame() %>%
        select(ends_with("_diff")) %>%
        colnames(.) %>%
        gsub("_diff", "", .)
      valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                  paste0(valid_cntrsts, collapse = "', '"),
                                  "'")
      stop("The contrast(s) is/are not valid, ",
           "please run `select_contrasts()` with a valid contrast as argument\n",
           valid_cntrsts_msg,
           call. = FALSE)
    }

    if(length(matches) == 1) {
      filtered <- object[row_data[,matches],]
    } else {
      filtered <- object[apply(row_data[,matches], 1, all)]
    }
  }
  return(filtered)
}


#' Generate a results table
#'
#' \code{get_results} generates a results table from a proteomics or RNA-seq dataset
#' on which differential analysis was performed.
#'
#' @param object SummarizedExperiment or DEGdata,
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
get_results <- function(object) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(object, "SummarizedExperiment"))

  row_data <- rowData(object, use.names = FALSE)
  # Show error if inputs do not contain required columns
  if(class(object) != "DEGdata" && (any(!c("name", "ID") %in% colnames(row_data)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(object)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(object)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }

  if(class(object) == "DEGdata"){
    row_data$name = rownames(object) %>% make.unique()
    if("SYMBOL" %in% colnames(row_data) & nrow(object@geneinfo)>1){
      row_data$ID = object@geneinfo[,1]
    }else{
      row_data$ID = row_data$name
    }
  }

  # Obtain average protein-centered enrichment values per condition
  row_data$mean <- rowMeans(assay(object), na.rm = TRUE)
  centered <- assay(object) - row_data$mean
  centered <- data.frame(centered) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(object)), by = "ID")
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

  # Join into a results table+
  # if(class(object) == "DEGdata"){
  #   ids <- ifelse("SYMBOL" %in% colnames(row_data) & nrow(object@geneinfo)>1,
  #                 data.frame(name = row_data$SYMBOL,
  #                            ID = object@geneinfo[,1]),
  #                 data.frame(name = rownames(row_data),
  #                            ID = rownames(row_data)))
  #
  # }else{
  if(check_ptm(object)){ # if is PTM with gene_name and protein_ID columns
    ids <- as.data.frame(row_data) %>% select(name, ID, gene_name, protein_ID, modified_aa, modified_pos)
  }else{ # pg data
    ids <- as.data.frame(row_data) %>% select(name, ID)
  }
  # }
  table <- left_join(ids, pval, by = c("name" = "rowname"))

  table <- left_join(table, df, by = c("name" = "rowname")) %>%
    arrange(desc(significant))

  return(table)
}

# Check whether a object is generated from a make_se_ptm and contain "gene_name", "protein_ID" colmuns.
check_ptm <- function(object){
  return(all(c("gene_name", "protein_ID") %in% colnames(rowData(object))))
}

# Internal function to obtain a table suitable for shiny visualization
get_table <- function(results, type = c("centered", "contrast")) {
  assertthat::assert_that(is.data.frame(results))
  type <- match.arg(type)

  if(length(grep("_ratio$", colnames(results))) < 1) {
    stop("'[contrast]_ratio' columns are not present in '",
         deparse(substitute(results)),
         "'.\nRun get_results() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_centered$", colnames(results))) < 1) {
    stop("'[contrast]_centered' columns are not present in '",
         deparse(substitute(results)),
         "'.\nRun get_results() to obtain the required columns",
         call. = FALSE)
  }

  # Only significant proteins
  significant <- results %>%
    filter(significant) %>%
    select(-significant)

  # Make centered table
  if(type == "centered") {
    cols <- grep("_ratio", colnames(significant))
    table <- significant[,-cols]
    colnames(table)[c(1,2)] <- c("Protein Name", "Protein ID")
    colnames(table)[grep("significant", colnames(table))] <-
      gsub("[.]", " - ", colnames(table)[grep("significant", colnames(table))])
    colnames(table) <- gsub("_centered", "", colnames(table)) %>% gsub("[_]", " ", .)
  }
  # Make contrast table
  if(type == "contrast") {
    cols <- grep("_centered", colnames(significant))
    table <- significant[,-cols]
    colnames(table)[c(1,2)] <- c("Protein Name", "Protein ID")
    colnames(table)[grep("significant", colnames(table))] <-
      gsub("[.]", " - ", colnames(table)[grep("significant", colnames(table))])
    colnames(table) <- gsub("_ratio", "", colnames(table)) %>% gsub("[_]", " ", .)
  }
  return(table)
}


# rm end as "\\.\\d*$" of genes used for combine panels
rm_digit_end <- function(x){
  gsub("\\.\\d*$", "", x)
}

