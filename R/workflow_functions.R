#' TMT workflow
#'
#' \code{TMT} is a wrapper function running the
#' entire differential enrichment/expression analysis workflow
#' for TMT-based proteomics data.
#' The protein table from
#' \href{http://www.nature.com/nprot/journal/v10/n10/full/nprot.2015.101.html}{IsobarQuant}
#' is used as direct input.
#'
#' @param proteins Data.frame,
#' The data object.
#' @param expdesign Data.frame,
#' The experimental design object.
#' @param fun "man", "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "min", "zero", "mixed" or "nbavg",
#' Function used for data imputation based on \code{\link{manual_impute}}
#' and \code{\link[MSnbase:impute-methods]{impute}}.
#' @param control Character(1),
#' The sample name to which the contrasts are generated
#' (the control sample would be most appropriate).
#' @param type 'all', 'control' or 'manual',
#' The type of contrasts that will be generated.
#' @param test Character,
#' The contrasts that will be tested if type = "manual".
#' These should be formatted as "SampleA_vs_SampleB" or
#' c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @param name Character(1),
#' Name of the column representing gene names.
#' @param ids 'Character(1),
#' Name of the column representing protein IDs.
#' @param alpha Numeric(1),
#' sets the false discovery rate threshold.
#' @param lfc Numeric(1),
#' sets the log fold change threshold.
#' @return A list of 8 objects:
#' \item{se}{SummarizedExperiment object containing the original data}
#' \item{filt}{SummarizedExperiment object containing the filtered data}
#' \item{norm}{SummarizedExperiment object containing the normalized data}
#' \item{imputed}{SummarizedExperiment object containing the imputed data}
#' \item{diff}{SummarizedExperiment object
#' containing FDR estimates of differential expression}
#' \item{dep}{SummarizedExperiment object
#' annotated with logical columns indicating significant proteins}
#' \item{results}{data.frame containing containing
#' all results variables from the performed analysis}
#' \item{param}{data.frame containing the test parameters}
#' @examples
#' \dontrun{
#'
#' TMT_res <- TMT()
#'
#' }
#' @export
TMT <- function(proteins, expdesign,
                fun = c("man", "bpca", "knn", "QRILC", "MLE", "MinDet",
                        "MinProb", "min", "zero", "mixed", "nbavg"),
                type = c("all", "control", "manual"), control = NULL, test = NULL,
                name = "gene_name", ids = "protein_id",
                alpha = 0.05, lfc = 1) {
    # Show error if inputs are not the required classes
    if(is.integer(alpha)) alpha <- as.numeric(alpha)
    if(is.integer(lfc)) lfc <- as.numeric(lfc)
    assertthat::assert_that(is.data.frame(proteins),
                            is.data.frame(expdesign),
                            is.character(fun),
                            is.character(type),
                            is.character(name),
                            length(name) == 1,
                            is.character(ids),
                            length(ids) == 1,
                            is.numeric(alpha),
                            length(alpha) == 1,
                            is.numeric(lfc),
                            length(lfc) == 1)

    # Show error if inputs do not contain required columns
    type <- match.arg(type)
    fun <- match.arg(fun)

    if(length(grep(paste("^", name, "$", sep = ""), colnames(proteins))) < 1) {
        stop("'", name, "' is not a column in '",
             deparse(substitute(proteins)), "'.",
             call. = FALSE)
    }
    if(length(grep(paste("^", ids, "$", sep = ""), colnames(proteins))) < 1) {
        stop("'", ids, "' is not a column in '",
             deparse(substitute(proteins)), "'.",
             call. = FALSE)
    }
    if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
        stop("'label', 'condition' and/or 'replicate' columns are ",
             "not present in the experimental design",
             call. = FALSE)
    }
    # Show error if inputs are not valid
    if(!type %in% c("all", "control")) {
        stop("run TMT() with a valid type",
             "\nValid types are: 'all', 'control' and 'manual'.",
             call. = FALSE)
    }

    # If input is a tibble, convert to data.frame
    if(tibble::is_tibble(proteins))
      proteins <- as.data.frame(proteins)
    if(tibble::is_tibble(expdesign))
      expdesign <- as.data.frame(expdesign)

    # Filter the proteins for Reverse hits (indicated by '###' in the gene_name)
    proteins <- proteins[-grep("###", proteins$gene_name), ]
    # Make unique names and turn the data into a SummarizedExperiment
    cols <- grep("signal_sum", colnames(proteins))  # The reporter signal columns
    proteins_unique <- make_unique(proteins, name, ids, delim = "[|]")
    se <- make_se(proteins_unique, cols, expdesign)
    # Filter on missing values
    filt <- filter_missval(proteins_unique)
    # Variance stabilization
    norm <- normalize_vsn(filt)
    # Impute missing values
    imputed <- impute(norm, fun)
    # Test for differential expression by empirical Bayes moderation
    # of a linear model and defined contrasts
    diff <- test_diff(imputed, type, control)
    # Denote differential expressed proteins
    dep <- add_rejections(diff, alpha, lfc)
    # Generate a results table
    results <- get_results(dep)

    param <- data.frame(alpha, lfc)
    return(list(se = se, filt = filt, norm = norm, imputed = imputed,
                diff = diff, dep = dep, results = results, param = param))
}

#' LFQ workflow
#'
#' \code{LFQ} is a wrapper function running the entire
#' differential enrichment/expression analysis workflow
#' for label free quantification (LFQ)-based proteomics data.
#' The protein table from
#' \href{http://www.nature.com/nbt/journal/v26/n12/full/nbt.1511.html}{MaxQuant}
#' is used as direct input.
#'
#' @param proteins Data.frame,
#' The data object.
#' @param expdesign Data.frame,
#' The experimental design object.
#' @param fun "man", "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "min", "zero", "mixed" or "nbavg",
#' Function used for data imputation based on \code{\link{manual_impute}}
#' and \code{\link[MSnbase:impute-methods]{impute}}.
#' @param control Character(1),
#' The sample name to which the contrasts are generated
#' (the control sample would be most appropriate).
#' @param type 'all', 'control' or 'manual',
#' The type of contrasts that will be generated.
#' @param test Character,
#' The contrasts that will be tested if type = "manual".
#' These should be formatted as "SampleA_vs_SampleB" or
#' c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @param filter Character,
#' Name(s) of the column(s) to be filtered on.
#' @param name Character(1),
#' Name of the column representing gene names.
#' @param ids 'Character(1),
#' Name of the column representing protein IDs.
#' @param alpha Numeric(1),
#' sets the false discovery rate threshold.
#' @param lfc Numeric(1),
#' sets the log fold change threshold.
#' @return A list of 9 objects:
#' \item{data}{data.frame containing the original data}
#' \item{se}{SummarizedExperiment object containing the original data}
#' \item{filt}{SummarizedExperiment object containing the filtered data}
#' \item{norm}{SummarizedExperiment object containing the normalized data}
#' \item{imputed}{SummarizedExperiment object containing the imputed data}
#' \item{diff}{SummarizedExperiment object
#' containing FDR estimates of differential expression}
#' \item{dep}{SummarizedExperiment object
#' annotated with logical columns indicating significant proteins}
#' \item{results}{data.frame containing containing
#' all results variables from the performed analysis}
#' \item{param}{data.frame containing the test parameters}
#' @examples
#'
#' data <- UbiLength
#' expdesign <- UbiLength_ExpDesign
#' results <- LFQ(data, expdesign, 'MinProb', 'control', 'Ctrl')
#'
#' @export
LFQ <- function(proteins, expdesign,
                fun = c("man", "bpca", "knn", "QRILC", "MLE", "MinDet",
                        "MinProb", "min", "zero", "mixed", "nbavg"),
                type = c("all", "control", "manual"), control = NULL, test = NULL,
                filter = c("Reverse", "Potential.contaminant"),
                name = "Gene.names", ids = "Protein.IDs",
                alpha = 0.05, lfc = 1) {
    # Show error if inputs are not the required classes
    if(is.integer(alpha)) alpha <- as.numeric(alpha)
    if(is.integer(lfc)) lfc <- as.numeric(lfc)
    assertthat::assert_that(is.data.frame(proteins),
                            is.data.frame(expdesign),
                            is.character(fun),
                            is.character(type),
                            is.character(filter),
                            is.character(name),
                            length(name) == 1,
                            is.character(ids),
                            length(ids) == 1,
                            is.numeric(alpha),
                            length(alpha) == 1,
                            is.numeric(lfc),
                            length(lfc) == 1)

    # Show error if inputs do not contain required columns
    if(length(grep(paste("^", name, "$", sep = ""), colnames(proteins))) < 1) {
        stop("'", name, "' is not a column in '",
             deparse(substitute(proteins)), "'.",
             call. = FALSE)
    }
    if(length(grep(paste("^", ids, "$", sep = ""), colnames(proteins))) < 1) {
        stop("'", ids, "' is not a column in '",
             deparse(substitute(proteins)), "'.",
             call. = FALSE)
    }
    if(length(grep(paste("^", filter, "$", sep = "", collapse = "|"),
                    colnames(proteins))) < 1) {
        stop("Not all filter columns are present in '",
             deparse(substitute(proteins)), "'.",
             call. = FALSE)
    }
    if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
        stop("'label', 'condition' and/or 'replicate' columns are ",
             "not present in the experimental design",
             call. = FALSE)
    }
    # Show error if inputs are not valid
    if (!type %in% c("all", "control")) {
        stop("run LFQ() with a valid type",
             "\nValid types are: 'all', 'control' and 'manual'.",
             call. = FALSE)
    }

    # If input is a tibble, convert to data.frame
    if(tibble::is_tibble(proteins)) proteins <- as.data.frame(proteins)
    if(tibble::is_tibble(expdesign)) expdesign <- as.data.frame(expdesign)

    # Filter out the positive proteins (indicated by '+')
    # in the pre-defined columns
    cols_filt <- grep(paste("^", filter, "$", sep = "", collapse = "|"),
                      colnames(proteins))  # The columns to filter on
    if (!is.null(cols_filt)) {
        if (length(cols_filt) == 1) {
          rows <- which(proteins[, cols_filt] == "+")
          if(length(rows) > 0) proteins <- proteins[-rows,]
        } else {
          rows <- which(apply(proteins[, cols_filt] == "+", 1, any))
          if(length(rows) > 0) proteins <- proteins[-rows,]
        }
    }

    # Make unique names and turn the proteins into a SummarizedExperiment
    cols <- grep("^LFQ", colnames(proteins))  # The LFQ intensity columns
    proteins_unique <- make_unique(proteins, name, ids, delim = ";")
    se <- make_se(proteins_unique, cols, expdesign)
    # Filter on missing values
    filt <- filter_missval(se)
    # Variance stabilization
    norm <- normalize_vsn(filt)
    # Impute missing values
    imputed <- impute(norm, fun)
    # Test for differential expression by empirical Bayes moderation
    # of a linear model and defined contrasts
    diff <- test_diff(imputed, type, control)
    # Denote differential expressed proteins
    dep <- add_rejections(diff, alpha, lfc)
    # Generate a results table
    results <- get_results(dep)

    param <- data.frame(alpha, lfc)
    return(list(data = proteins_unique, se = se, filt = filt, norm = norm,
                imputed = imputed, diff = diff, dep = dep,
                results = results, param = param))
}

#' Generate a markdown report
#'
#' \code{report} generates a report of the analysis performed
#' by \code{\link{TMT}} and \code{\link{LFQ}} wrapper functions.
#' Additionally, the results table is saved as a tab-delimited file.
#'
#' @param results List of SummarizedExperiment objects obtained
#' from the \code{\link{LFQ}} or \code{\link{TMT}} wrapper functions.
#' @return A \code{\link[rmarkdown:rmarkdown-package]{rmarkdown}} report is generated and saved.
#' Additionally, the results table is saved as a tab-delimited txt file.
#' @examples
#' \dontrun{
#'
#' data <- UbiLength
#' expdesign <- UbiLength_ExpDesign
#'
#' results <- LFQ(data, expdesign, 'MinProb', 'control', 'Ctrl')
#' report(results)
#'
#' }
#' @export
report <- function(results) {
    # Show error if inputs are not the required classes
    assertthat::assert_that(is.list(results))

    # Show error in case that the required objects
    # are not present in the list object 'results'
    if(any(!c("data", "se", "filt", "norm",
               "imputed", "diff", "dep",
               "results", "param") %in% names(results))) {
        stop("run report() with appropriate input generated by LFQ() or TMT()",
             call. = FALSE)
    }

    # Extract the objects used in the Rmarkdown report from the results object
    data <- results$se
    filt <- results$filt
    norm <- results$norm
    dep <- results$dep
    param <- results$param
    table <- results$results

    message("Render reports")
    # Render the Rmarkdown report
    wd <- paste(getwd(),
      "/Report_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      sep = "")
    dir.create(wd)
    file <- paste(system.file(package = "DEP2"), "/Report.Rmd", sep = "")
    rmarkdown::render(file, output_format = "all",
                      output_dir = wd, quiet = TRUE)

    message("Save tab-delimited table")
    # Save the results table in a tab-delimited txt file
    utils::write.table(table, paste(wd, "results.txt", sep = "/"),
                       row.names = FALSE, sep = "\t")

    message("Save RData object")
    # Save the results object for later use
    save(results, file = paste(wd, "results.RData", sep = "/"))

    message(paste0("Files saved in: '", wd, "'"))
}
