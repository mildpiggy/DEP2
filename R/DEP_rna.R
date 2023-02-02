
## creact a class to store DEseq2 result
#' Class "DEGdata"
#' This class store the test result from DESeq2
#'
#' @name DEGdata-class
#' @slot geneinfo data.frame. Store gene IDs and SYMBOL information.
#' @slot ntf matrix. Log2 transform assay
#' @slot rlg matrix. rlog transform assay
#' @slot test_result DataFrame. Test result from DESeq2
#' @importClassesFrom DESeq2 DESeqDataSet
#' @exportClass DEGdata
setClass("DEGdata",
         contains = "DESeqDataSet",
         slots = representation(
           geneinfo = "data.frame",
           ntf = "matrix",
           rlg = "matrix",
           test_result = "DataFrame")
)


setMethod("[", "DEGdata",
          function(x, i, j, drop = FALSE) {
            geneinfo <- x@geneinfo
            ntf <- x@ntf
            rlg <- x@rlg
            test_result <- x@test_result

            if (!missing(i)) {
              if (is.character(i)) {
                fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
                i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                  i, rownames(x), fmt
                )
              }
              i <- as.vector(i)
              if(nrow(geneinfo) == nrow(x)) geneinfo <- geneinfo[i,,drop=drop]
              if(nrow(ntf) == nrow(x)) ntf <- ntf[i,,drop=drop]
              if(nrow(rlg) == nrow(x)) rlg <- rlg[i,,drop=drop]
              if(nrow(test_result) == nrow(x)) test_result <- test_result[i,,drop=drop]
            }

            if (!missing(j)) {
              if (is.character(j)) {
                fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
                j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                  j, colnames(x), fmt
                )
              }
              j <- as.vector(j)
              if(ncol(ntf) == ncol(x)) ntf <- ntf[,j,drop=drop]
              if(ncol(rlg) == ncol(x)) rlg <- rlg[,j,drop=drop]
            }

            out <- callNextMethod()
            BiocGenerics:::replaceSlots(out, geneinfo=geneinfo, ntf=ntf,
                                        rlg=rlg, test_result=test_result,check=FALSE)
          }
)






# label = colnames(dds_1_save)
# get_exdesign_parse(label, mode = "char")

## Construct experiment design basic on colnames parse. Result table contain four columns: label,ID,condition,replicate
#' Construct experiment design basic on colnames parse
#'
#' @param label Characters, label names(column names) of assay
#' @param mode one of "delim", "char"
#' @param chars Integer(1). If mode is "char", how many chars in the end of label are 'replicate'
#' @param sep Character(1). If mode is "sep", how many chars in the end of label are 'replicate'
#'
#' @return
#'  A data.frame contain four columns: label,ID,condition,replicate
#' @export
#' @examples
#' ## Parse on delim
#' samples = c("Quantity.A_1","Quantity.A_2","Quantity.A_3","Quantity.B_1","Quantity.B_2","Quantity.B_3",
#'             "Quantity.C_1","Quantity.C_2")
#' get_exdesign_parse(samples, mode = "delim", sep = "_")
#'
#' ## Parse on certain character
#' samples = c("Quantity.A1","Quantity.A2","Quantity.A3","Quantity.B1","Quantity.B2","Quantity.B3",
#'             "Quantity.C1","Quantity.C2")
#' get_exdesign_parse(samples, mode = "char", chars = 1)
#' # Reserve prefix
#' get_exdesign_parse(samples, mode = "char", chars = 1, remove_prefix = F)
get_exdesign_parse <- function(label, mode = c("delim", "char"),
                               chars = 1, sep = "_", remove_prefix = T, remove_suffix = F){
  mode = match.arg(mode)

  if(remove_prefix){
    label <- delete_prefix(label) %>% make.names()
  }
  if(remove_suffix){
    label <- delete_suffix(label) %>% make.names()
  }

  if (mode == "char") {
    col_data <- data.frame(label = label, stringsAsFactors = FALSE) %>%
      mutate(condition = substr(label, 1, nchar(label) -
                                  chars), replicate = substr(label, nchar(label) +
                                                               1 - chars, nchar(label))) %>% unite(ID, condition,
                                                                                                   replicate, remove = FALSE)
  }
  if (mode == "delim") {
    # colnames(raw) = gsub(get_suffix(colnames(raw)),"", colnames(raw))
    col_data <- data.frame(label = label, stringsAsFactors = FALSE) %>%
      separate(label, c("condition", "replicate"), sep = sep,
               remove = FALSE, extra = "merge") %>% unite(ID,
                                                          condition, replicate, remove = FALSE)
  }
  # rownames(col_data) <- col_data$ID
  rownames(col_data) <- seq(nrow(col_data))
  return(col_data)
}



#' Differential expression test on a DESeqDataSet
#'
#' \code{test_diff_deg} performs a differential expression test based on
#' \code{DESeq2} package and format the test results from \link[DESeq2]{DESeq}() and \link[DESeq2]{results}().
#'
#' @param dds A DESeqDataSet object.
#'
#' @param type "control", "all" or "manual",
#' The type of contrasts that will be tested.
#' This can be all possible pairwise comparisons ("all"),
#' limited to the comparisons versus the control ("control"), or
#' manually defined contrasts ("manual").
#' @param control Character(1),
#' The condition to which contrasts are generated if type = "control"
#' (a control condition would be most appropriate).
#' @param test_contrasts Character,
#' The contrasts that will be tested if type = "manual".
#' These should be formatted as "SampleA_vs_SampleB" or
#' c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @param contrast_upon Character, the contrast to bulid upon which column in experimentdesign.
#' @param filter_ihw logical(1),Whether to use ihw to perform independent filtering and p-value adjustment, default FALSE. See \link[DESeq2]{results}
#' @param independentFiltering logical(1), whether independent filtering should be applied automatically, default TRUE. See \link[DESeq2]{results}
#' @param lfcshark logical(1),whether to use lfcshark, default FALSE. See \link[DESeq2]{results}
#' @param ... Parameters transmitted to \link[DESeq2]{DESeq2}
#' @export
#'
test_diff_deg <- function(dds, type = c("all", "control", "manual"), control = NULL,
                    test_contrasts = NULL, contrast_upon = "condition",
                    filter_ihw = FALSE, independentFiltering = TRUE, lfcshark = FALSE,
                    ... # DESeq parameters
                    ) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.null(control) | is.character(control),
                          is.null(test_contrasts) | is.character(test_contrasts),
                          is.character(contrast_upon)
  )
  type <- match.arg(type)

  coldata = colData(dds)
  conditions <- as.character(unique(coldata[ , contrast_upon]))

  cat("Performing DESeq analysis... \n")
  dds1 = DESeq2::DESeq(dds, ...)
  cat("Formatting analysis result... \n")
  dds2 = get_res(dds1, type = type, control = control, test_contrasts, contrast_upon = contrast_upon,
                 filter_ihw = filter_ihw, independentFiltering = independentFiltering, lfcshark = lfcshark)

  cat("Calulate a rlog & log2 transform assay, stored in rlg or ntd slot... \n")
  dds2 = rlg_deg(dds2)
  dds2 = ntf_deg(dds2)
  # # dds2@rlg <- rlog(dds1,blind=FALSE) %>% assay
  # ntd = normTransform(dds1)
  # if(!(all(rownames(assay(dds1)) == rownames(assay(ntd))) & all(rownames(assay(dds1)) == rownames(diff)))) {
  #   stop("the order is disrupted")
  # }
  # dds2@ntf = assay(ntd)

  rd <- rowData(dds2)
  rd = cbind(rd[,!(colnames(rd) %in% colnames(dds2@test_result))] , dds2@test_result)
  rowData(dds2) = rd
  return(dds2)
}



get_res <- function (dds, type = c("all", "control",  "manual"), control = NULL,
                     test_contrasts = NULL, contrast_upon = "condition", filter_ihw = FALSE,
                     independentFiltering = TRUE, lfcshark = FALSE)
{
  assertthat::assert_that(inherits(dds,"DESeqDataSet"),
                          is.character(type), length(type) ==1
                          )
  type <- match.arg(type)
  coldata = colData(dds)

  if(!contrast_upon %in% colnames(coldata))
    stop("contrast_upon:",contrast_upon,"is absent in experiement design.")

  conditions <- as.character(unique(coldata[, contrast_upon]))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                    collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control,
                                                    "- ", sep = " "), "", .) %>% paste(" - ", control,
                                                                                       sep = "")
      }
    }
  }
  if (type == "control") {
    assertthat::assert_that(is.character(control), length(control) == 1)
    if (is.null(control)){
      stop("Please select control")
    }
    if(!control %in% conditions)
      stop("control:", control,"do not exist in giving contrast_upon and dds.")
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = " - ")
  }
  if (type == "manual") {
    if (is.null(test_contrasts)) {
      stop("Please give one or more test_contrasts that you want to test, when type is 'manual'.")
    }
    assertthat::assert_that(is.character(test_contrasts))
    if (any(!unlist(strsplit(test_contrasts, "_vs_")) %in% conditions)) {
      stop("Valid contrasts should be like this: conditionA_vs_conditionB")
    }
    cntrst <- gsub("_vs_", " - ", test_contrasts)
  }
  message("extracted results: ", paste(gsub(" - ", "_vs_",
                                            cntrst), collapse = ", "))
  getres <- function(dds, contrast, filter_ihw = FALSE, independentFiltering = TRUE,
                     lfcshark = FALSE) {
    if (!filter_ihw) {
      res = results(dds, contrast = contrast, independentFiltering = independentFiltering)
    }
    else {
      res = results(dds, contrast = contrast, independentFiltering = independentFiltering,
                    filterFun = ihw)
    }
    if (lfcshark) {
      res = lfcShrink(dds = dds, contrast = contrast, res = res,
                      type = "ashr")
    }
    return(res)
  }
  res_df <- sapply(cntrst, function(i) {
    re <- getres(dds = dds, contrast = c(contrast_upon, gsub(" - .*",
                                                             "", i), gsub(".* - ", "", i)), filter_ihw = filter_ihw,
                 independentFiltering = independentFiltering, lfcshark = lfcshark)
  })
  names(res_df) = gsub(" - ", "_vs_", names(res_df))
  res_df = do.call(cbind, res_df)
  names(res_df) = gsub("\\.", "_", names(res_df))
  names(res_df) = gsub("_log2FoldChange", "_diff", names(res_df))
  names(res_df) = gsub("_pvalue", "_p.val", names(res_df))
  names(res_df) = gsub("_padj", "_p.adj", names(res_df))
  test_res <- new("DEGdata", dds)
  test_res@test_result = res_df
  return(test_res)
}

#' @rdname DEP2_add_rejections
#' @exportMethod add_rejections
setMethod("add_rejections",
          "DEGdata",
          function(diff,alpha,lfc,thresholdmethod,curvature,x0_fold){
            add_rejections.DEGdata(diff,alpha,lfc,thresholdmethod,curvature,x0_fold)
          }
)

add_rejections.DEGdata <- function (diff, alpha = 0.05, lfc = 1,thresholdmethod="intersect",curvature=1,x0_fold = 2)
{
  test_result <- diff@test_result
  # if (any(!c("name", "ID") %in% colnames(test_result))) {
  #   stop("'name' and/or 'ID' columns are not present in '",
  #        deparse(substitute(diff)), "'\nRun make_unique() and make_se() to obtain the required columns",
  #        call. = FALSE)
  # }
  if (length(grep("_p.adj|_diff", colnames(test_result))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(diff)), "'\nRun get_result() to obtain the required columns",
         call. = FALSE)
  }

  if( length(thresholdmethod)!=1 || (!thresholdmethod %in% c("intersect","curve")) ){
    stop("thresholdmethod should be one of intersect/curve")
  }
  if( length(grep("_significant|significant", colnames(test_result))) > 0 ){
    exit_sig_cols <- grep("_significant|significant", colnames(test_result))
    warning("Significant columns already exited in test_result.\n   Remove/replace values in columns ",
            paste(colnames(test_result)[exit_sig_cols],collapse ="; ") )
    test_result = test_result[,-exit_sig_cols]
  }

  cols_p <- grep("_p.adj", colnames(test_result))
  cols_diff <- grep("_diff", colnames(test_result))

  if(thresholdmethod == "intersect"){
    if (is.integer(alpha))
      alpha <- as.numeric(alpha)
    if (is.integer(lfc))
      lfc <- as.numeric(lfc)
    assertthat::assert_that(class(diff) == "DEGdata",
                            is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                            length(lfc) == 1)
    if (length(cols_p) == 1) {
      p_reject <- test_result[, cols_p] <= alpha
      p_reject[is.na(p_reject)] <- FALSE
      diff_reject <- abs(test_result[, cols_diff]) >= lfc
      diff_reject[is.na(diff_reject)] <- FALSE
      test_result$significant  <-
        test_result[,gsub("_p.adj", "_significant", colnames(test_result)[cols_p])] <-
        (p_reject & diff_reject)
      # test_result$contrast_significant <- test_result$significant
      # colnames(test_result)[ncol(test_result)] <- gsub("_p.adj", "_significant", colnames(test_result)[cols_p])
      diff@test_result = test_result
    }
    if (length(cols_p) > 1) {
      p_reject <- test_result[, cols_p] <= alpha
      p_reject[is.na(p_reject)] <- FALSE
      diff_reject <- abs(test_result[, cols_diff]) >= lfc
      diff_reject[is.na(diff_reject)] <- FALSE
      sign_df <- p_reject & diff_reject
      sign_df <- cbind(sign_df, significant = apply(sign_df,
                                                    1, function(x) any(x)))
      colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
      sign_df <- cbind(name = test_result$name, as.data.frame(sign_df))
      test_result <- left_join(test_result %>%as.data.frame(),
                               sign_df, by = "name")
      diff@test_result = test_result
    }
    rowData(diff) <- cbind(rowData(diff), diff@test_result)
    return(diff)
  }else if(thresholdmethod == "curve"){
    if (is.integer(curvature))
      curvature <- as.numeric(curvature)
    if (is.integer(x0_fold))
      x0_fold <- as.numeric(x0_fold)
    assertthat::assert_that(class(diff) == "DEGdata",
                            is.numeric(curvature), length(curvature) == 1, is.numeric(x0_fold),
                            length(x0_fold) == 1)

    if (length(cols_p) == 1) {
      x = test_result[, cols_diff]
      x = fun.outlier(x) %>% na.omit()
      fit <- fitnormal(x)
      Sigma = sqrt(fit$theta[2])
      x0 = x0_fold*Sigma
      polar = ifelse((test_result[,cols_diff] > 0),1,-1)
      contrast_significant = ( -log10(test_result[,cols_p]) > curvature/abs(test_result[,cols_diff] - polar*x0) &
                   ifelse((test_result[,cols_diff] >= 0),test_result[,cols_diff]>x0,test_result[,cols_diff] < -x0) )
      contrast_significant[is.na(contrast_significant)] = FALSE
      test_result$significant <-
        test_result[,gsub("_p.adj", "_significant", colnames(test_result)[cols_p])] <- contrast_significant

      # colnames(test_result)[ncol(test_result)] <- gsub("_p.adj",
      #                                                  "_significant", colnames(test_result)[cols_p])
      diff@test_result = test_result
    }
    if(length(cols_p) > 1){
      sign_df <- sapply(1:length(cols_p), function(i) {
        x = test_result[, cols_diff[i]]
        x = fun.outlier(x) %>% na.omit()
        fit <- fitnormal(x)
        Sigma = sqrt(fit$theta[2])
        x0 = x0_fold*Sigma
        polar = ifelse((test_result[,cols_diff[i]] > 0),1,-1)
        contrast_significant <- ( -log10(test_result[,cols_p[i]]) > curvature/abs(test_result[,cols_diff[i]] - polar*x0) &
                                    ifelse((test_result[,cols_diff[i]] >= 0),test_result[,cols_diff[i]]>x0,test_result[,cols_diff[i]] < -x0) )
        contrast_significant[is.na(contrast_significant)] =FALSE
      })
      colnames(sign_df) =  gsub("_p.adj", "significant", colnames(test_result)[cols_p])
      sign_df <- cbind(sign_df, significant = apply(sign_df,
                                                    1, function(x) any(x)))
      sign_df <- cbind(name = test_result$name, as.data.frame(sign_df))
      test_result <- left_join(test_result %>% as.data.frame(),
                               sign_df, by = "name")
      diff@test_result = test_result
      diff@test_result = test_result
    }

    rowData(diff) <- cbind(diff[,!(colnames(diff) %in% colnames(diff@test_result))] , diff@test_result)
    return(diff)
  }else stop("thresholdmethod should be one of 'intersect'/'curve'",
             call. = FALSE)

}


#' Assign normalized assay for DEGdata object
#'
#' \code{ntf_deg} transform a log2(default) transformed normalized assay based on \code{\link[DESeq2]{normTransform}()}
#'
#' @param DEGdata a DEGdata object from get_res
#' @inheritParams DESeq2::normTransform
#'
#' @return
#' A DEGdata with a normalized assay in ntf slot
#' @export
#'
#' @importFrom DESeq2 normTransform
ntf_deg <- function(DEGdata, f = log2, pc = 1){
  assertthat::assert_that(inherits(DEGdata,"DEGdata"))
  DEGdata@ntf = assay(normTransform(DEGdata, f = f, pc = pc))
  return(DEGdata)
}


#' Assign a rlog transformed assay for DEGdata object
#'
#' \code{rlg_deg} calulate a rolg transformed assay based on \code{\link[DESeq2]{rlog}()}
#'
#' @param DEGdata a DEGdata object from get_res
#' @inheritParams DESeq2::rlog
#' @param ... Other options in \code{\link[DESeq2]{rlog}()}.
#'
#' @return
#' A DEGdata with a rlog transform assay in rlg slot
#' @export
#'
rlg_deg <- function(DEGdata, blind = FALSE, ...){
  DEGdata@rlg <- rlog(DEGdata,blind, ...) %>% assay
  return(DEGdata)
}


## RNAseq workflow
'
count_mat <- read.csv("example/RNAseq_example_countmatrix3.txt", sep = "\t", stringsAsFactors = F, header = T)
count_mat <- count_mat %>% `rownames<-`(.[,1]) %>% .[,-1]
dds_1 <- DESeqDataSetFromMatrix(countData = count_mat,colData = get_exdesign_parse(colnames(count_mat)),
                                design = as.formula(paste0("~", paste(c("condition","replicate"), collapse = " + "))))

dds_filt <- DEP2::filter_se(dds_1, missnum = 2, thr = 1, rowsum_threshold = 25)
dds3 <- test_diff_deg(dds_filt)
dds <- DESeq(dds_filt)
colData(dds)
results(dds) %>% class
dds_res <- get_res(dds = dds,
                   type = "control",
                   control = "untrt",
                   test = NULL,
                   contrast_upon = "condition")
dds_res <- DEP2::add_rejections(dds_res)
dds_res = dds_res[complete.cases(dds_res@test_result),]
dds_res = rlg_deg(dds_res)
dds_res = ntf_deg(dds_res)




dds_res@test_result$trt_vs_untrt_diff %>% summary
plot_heatmap(normTransform(dds_res), label_few_peptide_rows = T)
plot_heatmap.DEGdata(dds_res)
plot_volcano(normTransform(dds_res) , contrast = get_contrast(dds_res)[1])
dds_res2 = ID_transform(dds_res)
dds_res2@geneinfo
'

