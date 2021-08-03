#' Import from MaxQuant
#'
#' \code{import_MaxQuant} imports a protein table from MaxQuant
#' and converts it into a SummarizedExperiment object.
#'
#' @param proteins Data.frame,
#' Protein table originating from MaxQuant.
#' @param expdesign Data.frame,
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' See \code{\link{UbiLength_ExpDesign}} for an example experimental design.
#' @param filter Character,
#' Name of the column(s) containing features to be filtered on.
#' @param intensities Character(1),
#' Prefix of the columns containing sample intensities.
#' @param names Character(1),
#' Name of the column containing feature names.
#' @param ids Character(1),
#' Name of the column containing feature IDs.
#' @param delim Character(1),
#' Sets the delimiter separating the feature names within on protein group.
#' @return A SummarizedExperiment object with log2-transformed values and
#' "name" and "ID" columns containing unique names and identifiers.
#' @examples
#' # Load example data and experimental design
#' data <- UbiLength
#' exp_design <- UbiLength_ExpDesign
#'
#' # Import data
#' se <- import_MaxQuant(data, exp_design)
#' @export
import_MaxQuant <- function(proteins,
                            expdesign,
                            filter = c("Reverse", "Potential.contaminant"),
                            intensities = "LFQ",
                            names = "Gene.names",
                            ids = "Protein.IDs",
                            delim = ";") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins),
                          is.data.frame(expdesign),
                          is.character(filter),
                          is.character(intensities),
                          length(intensities) == 1,
                          is.character(names),
                          length(names) == 1,
                          is.character(ids),
                          length(ids) == 1,
                          is.character(delim),
                          length(delim) == 1)

  # Show error if inputs do not contain required columns
  if(any(!filter %in% colnames(proteins))) {
    stop(paste0("'",
                paste0(filter, collapse = "' and/or '")
                , "' not found in '",
                deparse(substitute(proteins)), "'."),
         call. = FALSE)
  }

  columns <- grep(paste0("^", intensities), colnames(proteins))
  if(length(columns) < 2) {
    stop("specified intensities prefix ('", intensities,
         "') does not indicate >1 columns",
         "\nRun import_MaxQuant() with the appropriate intensities argument",
         call. = FALSE)
  }
  if(any(!apply(proteins[, columns], 2, is.numeric))) {
    stop("specified intensities prefix ('", intensities,
         "') does not indicate numeric columns",
         "\nRun import_MaxQuant() with the appropriate intensities argument",
         call. = FALSE)
  }

  # Filter based on 'filter' columns
  proteins <- filter_MaxQuant(proteins, filter)

  # Make unique names
  message("Making unique names")
  proteins_unique <- make_unique(proteins, names, ids, delim)

  # Obtain SE object
  message("Obtaining SummarizedExperiment object")
  make_se(proteins_unique, columns, expdesign)
}

#' Import from IsobarQuant
#'
#' \code{import_IsobarQuant} imports a protein table from IsobarQuant
#' and converts it into a SummarizedExperiment object.
#'
#' @param proteins Data.frame,
#' Protein table for which unique names will be created.
#' @param expdesign Data.frame,
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' See \code{\link{UbiLength_ExpDesign}} for an example experimental design.
#' @param intensities Character(1),
#' Prefix of the columns containing sample intensities.
#' @param names Character(1),
#' Name of the column containing feature names.
#' @param ids Character(1),
#' Name of the column containing feature IDs.
#' @param delim Character(1),
#' Sets the delimiter separating the feature names within on protein group.
#' @return A SummarizedExperiment object with log2-transformed values and
#' "name" and "ID" columns containing unique names and identifiers.
#' @examples
#' \dontrun{
#' # Load data
#' isobarquant_table <- read.csv("testfile.txt", header = TRUE,
#'                               stringsAsFactors = FALSE, sep = "\t")
#' exp_design <- read.csv("test_experimental_design.txt", header = TRUE,
#'                               stringsAsFactors = FALSE, sep = "\t")
#' # Import data
#' se <- import_IsobarQuant(isabarquant_table, exp_design)
#'
#' }
#' @export
import_IsobarQuant <- function(proteins,
                               expdesign,
                               intensities = "signal_sum",
                               names = "gene_name",
                               ids = "protein_id",
                               delim = "[|]") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins),
                          is.data.frame(expdesign),
                          is.character(intensities),
                          length(intensities) == 1,
                          is.character(names),
                          length(names) == 1,
                          is.character(ids),
                          length(ids) == 1,
                          is.character(delim),
                          length(delim) == 1)

  columns <- grep(paste0("^", intensities), colnames(proteins))
  if(length(columns) < 2) {
    stop("specified intensities prefix ('", intensities,
         "') does not indicate >1 columns",
         "\nRun import_IsobarQuant() with the appropriate intensities argument",
         call. = FALSE)
  }
  if(any(!apply(proteins[, columns], 2, is.numeric))) {
    stop("specified intensities prefix ('", intensities,
         "') does not indicate numeric columns",
         "\nRun import_IsobarQuant() with the appropriate intensities argument",
         call. = FALSE)
  }

  # Filter decoy hits
  message("Filtering decoy hits")
  rows <- grep("###", proteins[, grep(paste0("^", names, "$"), colnames(proteins))])
  proteins <- proteins[-rows, ]

  # Make unique names
  message("Making unique names")
  proteins_unique <- make_unique(proteins, names, ids, delim)

  # Obtain SE object
  message("Obtaining SummarizedExperiment object")
  make_se(proteins_unique, columns, expdesign)
}

#' Proteomics data processing
#'
#' \code{process} performs data processing on a SummarizedExperiment object.
#' It (1) filters a proteomics dataset based on missing values,
#' (2) applies variance stabilizing normalization and
#' (3) imputes eventual remaining missing values.
#'
#' @param se SummarizedExperiment,
#' Proteomics data with unique names and identifiers
#' annotated in 'name' and 'ID' columns.
#' The appropriate columns and objects can be generated
#' using the wrapper import functions \code{\link{import_MaxQuant}}
#' and \code{\link{import_IsobarQuant}} or the generic functions
#' \code{\link{make_se}} and \code{\link{make_se_parse}}.
#' @param thr Integer(1),
#' Sets the threshold for the allowed number of missing values per condition.
#' @param fun "man", "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "min", "zero", "mixed" or "nbavg",
#' Function used for data imputation based on \code{\link{manual_impute}}
#' and \code{\link[MSnbase:impute-methods]{impute}}.
#' @param ... Additional arguments for imputation functions as depicted in
#' \code{\link{manual_impute}} and \code{\link[MSnbase:impute-methods]{impute}}.
#' @return A filtered, normalized and imputed SummarizedExperiment object.
#' @examples
#' # Load datasets
#' data <- UbiLength
#' exp_design <- UbiLength_ExpDesign
#'
#' # Import data
#' se <- import_MaxQuant(data, exp_design)
#'
#' # Process data
#' processed <- process(se)
#' @export
process <- function(se, thr = 0, fun = c("man", "bpca", "knn", "QRILC", "MLE",
                                         "MinDet", "MinProb", "min", "zero",
                                         "mixed", "nbavg"), ...) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(thr),
                          length(thr) == 1,
                          is.character(fun))
  fun <- match.arg(fun)

  # Filter on missing values
  message("Filtering on missing values")
  filtered <- filter_missval(se, thr)

  # Normalize using vsn
  message("Normalizing data")
  normalized <- normalize_vsn(filtered)
  p1 <- plot_normalization(filtered, normalized)

  # Impute missing values
  message("Imputing remaining missing values")
  imputed <- impute(normalized, fun, ...)

  print(p1)
  return(imputed)
}

#' Differential expression analysis
#'
#' \code{analyze_dep} tests for differential expression of proteins
#' based on protein-wise linear models and empirical Bayes
#' statistics using \pkg{limma}.
#'
#' @param se SummarizedExperiment,
#' Proteomics data with unique names and identifiers
#' annotated in 'name' and 'ID' columns.
#' Additionally, the colData should contain sample annotation including
#' 'label', 'condition' and 'replicate' columns.
#' The appropriate columns and objects can be generated
#' using \code{\link{make_se}} or \code{\link{make_se_parse}}.
#' @param control Character(1),
#' The condition to which contrasts are generated
#' (a control condition would be most appropriate).
#' @param type "all", "control" or "manual",
#' The type of contrasts that will be tested.
#' This can be all possible pairwise comparisons ("all"),
#' limited to the comparisons versus the control ("control"), or
#' manually defined contrasts ("manual").
#' @param alpha Numeric(1),
#' Sets the threshold for the adjusted P value.
#' @param lfc Numeric(1),
#' Sets the threshold for the log2 fold change.
#' @param test Character,
#' The contrasts that will be tested if type = "manual".
#' These should be formatted as "SampleA_vs_SampleB" or
#' c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @param design_formula Formula,
#' Used to create the design matrix.
#' @return A SummarizedExperiment object
#' containing FDR estimates of differential expression and
#' logical columns indicating significant proteins.
#' @examples
#' # Load datasets
#' data <- UbiLength
#' exp_design <- UbiLength_ExpDesign
#'
#' # Import and process data
#' se <- import_MaxQuant(data, exp_design)
#' processed <- process(se)
#'
#' # Differential protein expression analysis
#' dep <- analyze_dep(processed, "control", "Ctrl")
#' dep <- analyze_dep(processed, "control", "Ctrl",
#'     alpha = 0.01, lfc = log2(1.5))
#' dep <- analyze_dep(processed, "manual", test = c("Ubi6_vs_Ubi4"))
#' @export
analyze_dep <- function(se, type = c("all", "control", "manual"),
                     control = NULL, alpha = 0.05, lfc = 1,
                     test = NULL, design_formula = formula(~ 0 + condition)) {

  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          is.numeric(alpha),
                          length(alpha) == 1,
                          is.numeric(lfc),
                          length(lfc) == 1,
                          class(design_formula) == "formula")
  type <- match.arg(type)

  # Test for differentially enriched proteins
  message("Test for differentially enriched proteins")
  diff <- test_diff(se, type, control, test, design_formula)

  # Add rejections
  message("Add rejections")
  add_rejections(diff, alpha, lfc)
}

#' Visualize the results in different types of plots
#'
#' \code{plot_all} visualizes the results of the differential protein
#' expression analysis in different types of plots. These are (1)
#' volcano plots, (2) heatmaps, (3) single protein plots, (4) frequency plots
#' and/or (5) comparison plots.
#'
#' @param dep SummarizedExperiment,
#' Data object which has been generated by \code{\link{analyze_dep}}
#' or the combination of \code{\link{test_diff}}
#' and \code{\link{add_rejections}}.
#' @param plots "volcano", "heatmap", "single", "freq" and/or "comparison",
#'
#' @return Pdfs containg the desired plots.
#' @examples
#' # Load datasets
#' data <- UbiLength
#' exp_design <- UbiLength_ExpDesign
#'
#' # Import and process data
#' se <- import_MaxQuant(data, exp_design)
#' processed <- process(se)
#'
#' # Differential protein expression analysis
#' dep <- analyze_dep(processed, "control", "Ctrl")
#'
#' \dontrun{
#' # Plot all plots
#' plot_all(dep)
#' }
#' @export
plot_all <- function(dep, plots = c("volcano", "heatmap",
                                    "single", "freq", "comparison")) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(plots))

  # Get row data
  row_data <- rowData(dep, use.names = FALSE)

  # Show error if inputs do not contain required columns
  dep_mcolnames <- colnames(row_data)
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun make_unique() to obtain required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }
  possible_plots <- c("volcano", "heatmap",
                      "single", "freq", "comparison")
  if(any(!plots %in% possible_plots)) {
    stop("Run plot_all() with a valid plots.",
         deparse(substitute(dep)),
         "'\nValid plots are one or combinations of: '",
         paste0(possible_plots, collapse = "', '"), "'.",
         call. = FALSE)
  }

  # Get all contrasts
  contrasts <- row_data %>%
    data.frame() %>%
    select(ends_with("_diff")) %>%
    colnames(.) %>%
    gsub("_diff", "", .)

  # Get all differentially enriched proteins
  significants <- row_data %>%
    data.frame() %>%
    filter(significant) %>%
    select(name)

  # Number of significants
  num <- nrow(significants)


  if("volcano" %in% plots) {
    # Plot volcanos
    message("Plot volcano:")
    grDevices::pdf(file = "plot_volcanos.pdf")
    for(contrast in contrasts) {
      message(contrast)
      print(plot_volcano(dep, contrast) + labs(title = contrast))
    }
    grDevices::dev.off()
  }

  if("heatmap" %in% plots) {
    # Plot heatmaps
    message("Plot heatmap:")
    grDevices::pdf(file = "plot_heatmaps.pdf", height = 10)
    message("Centered")
    plot_heatmap(dep, type = "centered",
                 row_font_size = ifelse(num > 400, 1, 400/num))
    message("Contrast")
    plot_heatmap(dep, type = "contrast",
                 row_font_size = ifelse(num > 400, 1, 400/num))
    grDevices::dev.off()
  }

  if("single" %in% plots) {
    # Plot individual plots
    message("Plot single proteins:")
    grDevices::pdf(file = "plot_single_proteins.pdf",
                   height = 10, width = 5)
    for(name in significants$name) {
      message(name)
      p1 <- plot_single(dep, name, type = "centered")
      p2 <- plot_single(dep, name, type = "contrast")
      gridExtra::grid.arrange(p1, p2, ncol = 1)
    }
    grDevices::dev.off()
  }

  if("freq" %in% plots) {
    # Plot frequency graphs
    message("Plot frequency graphs")
    grDevices::pdf(file = "plot_frequency_graphs.pdf")
    print(plot_cond_freq(dep))
    plot_cond_overlap(dep)
    plot_cond(dep)
    grDevices::dev.off()
  }

  if("comparison" %in% plots) {
    # Plot comparison graphs
    message("Plot comparison graphs")
    grDevices::pdf(file = "plot_comparison_graphs.pdf")
    print(plot_pca(dep, n = ifelse(nrow(dep) > 500, 500, nrow(dep))))
    plot_cor(dep)
    grDevices::dev.off()
  }

}

# Internal function to filter on specified columns
filter_MaxQuant <- function(proteins, filter_column_names) {
  assertthat::assert_that(is.data.frame(proteins),
                          is.character(filter_column_names))

  # Get columns
  cols_filt <- match(filter_column_names, colnames(proteins))

  # Check columns
  if(all(is.na(cols_filt))) {
    warning("No filtering applied\nSpecified filter_column_names ('",
            paste0(filter_column_names, collapse = "' and '"),
            "') do not indicate any column",
            call. = FALSE)
  }
  if(any(is.na(cols_filt))) {
    cols_filt <- cols_filt[!is.na(cols_filt)]
  }

  # Filter proteins based on 'filter_column_names' columns
  message("Filtering based on '", paste(filter_column_names, collapse = "', '"), "' column(s)")
  if (!is.null(cols_filt)) {
    NAs <- is.na(proteins[,cols_filt])
    proteins[,cols_filt][NAs] <- ""
    if (length(cols_filt) == 1) {
      proteins <- filter(proteins, proteins[,cols_filt] != "+")
    } else if(length(cols_filt) > 1) {
      proteins <- filter(proteins, !apply(proteins[,cols_filt] == "+", 1, any))
    }
  }
  return(proteins)
}

# Internal function to exclude differentially expressed proteins from
# particular contrasts
exclude_deps <- function(dep, contrasts) {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))

  row_data <- rowData(dep, use.names = FALSE)
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }

  if(is.null(contrasts)) {
    filtered <- dep
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
           "please run `exclude_deps()` with a valid contrast as argument\n",
           valid_cntrsts_msg,
           call. = FALSE)
    }

    if(length(matches) == 1) {
      filtered <- dep[!row_data[,matches],]
    } else {
      filtered <- dep[!apply(row_data[,matches], 1, any)]
    }
  }
  return(filtered)
}

# Internal function to select differentially expressed proteins from
# particular contrasts
select_deps <- function(dep, contrasts) {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))

  row_data <- rowData(dep, use.names = FALSE)
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun add_rejections() to obtain the required columns",
         call. = FALSE)
  }

  if(is.null(contrasts)) {
    filtered <- dep
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
           "please run `exclude_deps()` with a valid contrast as argument\n",
           valid_cntrsts_msg,
           call. = FALSE)
    }

    if(length(matches) == 1) {
      filtered <- dep[row_data[,matches],]
    } else {
      filtered <- dep[apply(row_data[,matches], 1, all)]
    }
  }
  return(filtered)
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

