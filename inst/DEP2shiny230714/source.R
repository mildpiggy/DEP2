library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(limma)
library(colourpicker)
library(tibble)
library(RColorBrewer)
library(BiocParallel)
library(DESeq2)
library(dplyr)


# Normalize SE object, like normalize_pe.
normalize_se = function(se,
                        method = c("vsn", "diff.median", "quantiles", "quantiles.robust")
){
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  method = match.arg(method)

  if(method == "vsn"){
    norm = normalize_vsn(se)
  }else{
    norm = se
    mt = assay(norm)
    mt = normalize_matrix(mt, method)
    assay(norm) = mt
  }
  return(norm)

}

# Adjuct protein intensity though standard proteins, like spike-in proteins.
correct_AP_by_protein = function(object, standard_protein,
                                 correct_key = c("name", "ID")
){
  assertthat::assert_that(inherits(object, "SummarizedExperiment"),
                          is.character(standard_protein), is.character(correct_key))

  mt = assay(object)
  rd = rowData(object)

  correct_key = match.arg(correct_key)

  standard_protein = standard_protein[standard_protein %in% rd[,correct_key]]
  if(length(standard_protein)<1){
    stop("The standart_protein do not exist in the object! Please check your input!")
  }
  message("The assay is corrected based on proteins ",paste0(standard_protein,collapse = ", "))

  mt_stander = mt[which(rd[,correct_key] %in% standard_protein),]
  standarmean = rowMeans(mt_stander)
  correct_coefficient = colMeans(mt_stander - standarmean)

  mt_corrected = t(t(mt) - correct_coefficient)

  assay(object) = mt_corrected
  return(object)
}

# Clean expression colnames, make sure labels are suitable for get_exdesign_parse
clean_ecolnames = function(label){
  endpattern = paste0("\\d+",get_suffix(label),"+$")

  # remove '-' by "_"
  label = gsub("-", "-", label, perl=TRUE)
  # if it is '.' that connect replication, replacement '.' by' _'. like 'A.1' to 'A_2'
  label = gsub(paste0("\\.(?=",endpattern,")"), "_", label, perl=TRUE)
  # remove '_' except the last one which link last numbers (replication)
  label = gsub(paste0("\\_(?!",endpattern,")"), ".", label, perl=TRUE)
  # remove '_' except the last one which link last numbers (replication)
  label = gsub(paste0("_(?!",endpattern,")"), ".", label, perl=TRUE)

  return(label)
}




## for RNAseq differential analysis



#' Plot normalized counts for a gene
#'
#' Plot for normalized counts of a single gene, with jittered points superimposed
#' on the boxplot
#'
#' Note: this function relies on the \code{\link{plotCounts}} function of DESeq2,
#' therefore pseudocounts of 0.5 are added to each point
#'
#' @param dds A \code{\link{DESeqDataSet}} object.
#' @param gene A character, specifying the name of the gene to plot
#' @param intgroup Interesting groups: a character vector of
#' names in \code{colData(dds)} to use for grouping
#' @param annotation_obj A \code{data.frame} object, with \code{row.names} as gene
#' identifiers (e.g. ENSEMBL ids) and a column, \code{gene_name}, containing
#' e.g. HGNC-based gene symbols. Optional.
#' @param transform Logical value, corresponding whether to have log scale y-axis
#' or not. Defaults to TRUE.
#' @param ylimZero Logical value, if set the lower ylim to 0
#' @param add_labels if add_labels
#' @param labels_repel Logical value. Whether to use \code{ggrepel}'s functions to
#' place labels; defaults to TRUE.act when add_labels is TRUE
#'
#' @return An object created by \code{ggplot}
#' @export
ggplot_Counts <- function(dds, gene, intgroup = "condition", annotation_obj = NULL,
                          transform = TRUE, ylimZero = TRUE, add_labels = F, labels_repel = T) {
  df <- plotCounts(dds, gene, intgroup, returnData = TRUE)
  df$sampleID <- rownames(df)

  if (!is.null(annotation_obj)) {
    genesymbol <- annotation_obj$symbol[match(gene, rownames(annotation_obj))]
  } else {
    genesymbol <- ""
  }

  onlyfactors <- df[, match(intgroup, colnames(df))]
  df$plotby <- interaction(onlyfactors)
  # rownames(df) =

  p <-
    ggplot(df, aes_string(x = "plotby", y = "count", col = "plotby")) +#
    geom_boxplot(outlier.shape = NA, fill = "grey", colour = "black") +
    # geom_text(data = jittered_df,aes(x=conditionj,y=countj,label=sampleID)) +
    scale_x_discrete(name = "") +
    geom_jitter(aes_string(x = "plotby", y = "count"),
                position = position_jitter(width = 0.1), shape = 18, size = 5
    ) +
    scale_color_discrete(name = "Experimental\nconditions")

  if(add_labels){
    if (labels_repel) {
      p <- p + ggrepel::geom_text_repel(aes_string(label = "sampleID"), size = 4, col = "black")
    } else {
      p <- p + geom_text(aes_string(label = "sampleID"), hjust = -.1, vjust = 0, size = 4, col = "black")
    }
  }

  if (transform) {
    p <- p + scale_y_log10(name = "Normalized counts (log10 scale)")
  } else {
    p <- p + scale_y_continuous(name = "Normalized counts")
  }

  p <- p + theme_bw()

  if (!is.null(annotation_obj)) {
    p <- p + labs(title = paste0("Normalized counts for ", genesymbol, " - ", gene))
  } else {
    p <- p + labs(title = paste0("Normalized counts for ", gene))
  }

  if(ylimZero) {
    p <- p + ylim(0.1, NA)
  }

  p
}


# # #### for PR-Heatmap Connected to the whole app
# # # /**
# # #  * [description]
# # #  * @param  {[type]} rna [the full_results for downloading from DEG-RNAseq panel]
# # #  * @return {[type]}     [description]
# # #  */
# # prepare_rna <- function(rna) {
# #   names(rna)[1] = "name"
# #   names(rna) = paste("rna", names(rna), sep = "_")
# #   return(rna)
# # }
#
# # /**
# #  * [description]
# #  * @param  {[type]} pro [the full_dataset for downloading from DEP-LFQ panel]
# #  * @return {[type]}     [description]
# #  */
# prepare_pro <- function(pro) {
#   index_1 = which(gsub("_[[:digit:]]*", "", names(pro)) %in% unique(unlist(strsplit(gsub("_significant", "", grep("_significant", names(pro), value = T)), split = "_vs_"))))
#   index_2 = grep("_vs_|significant", names(pro))
#   pro = pro[ , c(1, index_1, index_2)]
#   names(pro) = paste("pro", names(pro), sep = "_")
#   return(pro)
# }




############# PPI function #############
library(data.table)
# library(visNetwork)
library(igraph)
library(htmlwidgets)
get_string_Env <- function () {
  if (!exists(".string_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".string_Env", new.env(), envir=envir)
  }
  get(".string_Env", envir = .GlobalEnv)
}




############# network function for PPI #############
##
stringNetwork <- function(linksTable,layoutway,nodecolor ,nodeshape,linecolor,nodesize,changesize=FALSE,
                          fontsize=200,changewidth=FALSE,linewidth,smoothline=FALSE,smoothtype="continous",highlightkey=TRUE){
  links3 <- linksTable
  nodes <- linksTable %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
  nodes$shape = nodeshape
  net <- igraph::graph_from_data_frame(d=links3,vertices=nodes,directed = F)
  if(changewidth){
    links3$width <- links3$combined_score/1500*linewidth
  }else{
    links3$width =linewidth
  }
  if(changesize){
    nodes$size <- igraph::degree(net)*nodesize/2
  }else{
    nodes$size =nodesize
  }
  if(highlightkey){
    nodes$color <- ifelse(igraph::degree(net)>=6,"#B54434",nodecolor)
  }else{
    nodes$color = nodecolor
  }
  # nodes$size =200
  colnames(nodes)[1] <- "id"

  nwplot <- visNetwork::visNetwork(nodes,links3,width = "170%",height = "500px") %>%
    visIgraphLayout(layout = layoutway ) %>%
    visNodes(size=nodes$size,
             # color = nodecolor ,
             font= list(size= fontsize,align = "left"),
             shadow = list(enabled = TRUE, size = 10)) %>%
    visEdges(shadow = F,
             color = list(color = linecolor, highlight = "red"),
             smooth = list(enabled=smoothline,type=smoothtype)
    ) %>%
    visOptions(highlightNearest = list(enabled = T, hover = T),
               nodesIdSelection = T)%>%
    visPhysics(solver = "forceAtlas2Based",
               forceAtlas2Based = list(gravitationalConstant = -500))

  return(nwplot)
}


# rm end as "\\.\\d*$" of genes used for combine panels
rm_digit_end <- function(x){
  gsub("\\.\\d*$", "", x)
}


################## for gene list tool

get_dep_genelist_Env <- function () {
  if (!exists(".genelist_dep_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_dep_Env", new.env(), envir=envir)
  }
  get(".genelist_dep_Env", envir = .GlobalEnv)
}

get_deg_genelist_Env <- function () {
  if (!exists(".genelist_deg_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_deg_Env", new.env(), envir=envir)
  }
  get(".genelist_deg_Env", envir = .GlobalEnv)
}

get_imported_genelist_Env <- function () {
  if (!exists(".genelist_imported_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_imported_Env", new.env(), envir=envir)
  }
  get(".genelist_imported_Env", envir = .GlobalEnv)
}

get_all_genelist_Env <- function () {
  if (!exists(".genelist_all_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_all_Env", new.env(), envir=envir)
  }
  get(".genelist_all_Env", envir = .GlobalEnv)
}

# assign DEP-LFQ and DEG-RNAseq all contrasts with log2 fc to this envir in order to import to GSEA options
get_gsea_genelist_Env <- function () {
  if (!exists(".genelist_gsea_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".genelist_gsea_Env", new.env(), envir=envir)
  }
  get(".genelist_gsea_Env", envir = .GlobalEnv)
}



statistics_plot <- function(dep,x = c("-log10.padj", "-log10.pval","padj","pval","t") ,y = c("-log10.padj", "-log10.pval","padj","pval","t"), contrast, plottype = c("x-y","histogram") ){
  row_data <- rowData(dep, use.names = FALSE)

  row_data = as.data.frame(row_data)
  p.adjcol <- grep(paste("^",contrast, "_p.adj", sep = ""),colnames(row_data))
  p.valcol <- grep(paste("^",contrast, "_p.val", sep = ""),colnames(row_data))
  tcol <- grep(paste("^",contrast, "_t.stastic", sep = ""),colnames(row_data))
  sigcol <- grep(paste(contrast, "_significant", sep = ""),colnames(row_data))
  df = data.frame(log10.padj = -log10(row_data[,p.adjcol]), log10.pval = -log10(row_data[,p.valcol]),
                  padj = row_data[,p.adjcol], pval = row_data[,p.valcol],
                  t = row_data[,tcol], significant = row_data[,sigcol], assay(dep))
  # if(x == "-log10.padj") {x = "log10.padj"}
  # if(x == "-log10.pval") {x = "log10.pval"}
  # if(y == "-log10.padj") {y = "log10.padj"}
  # if(y == "-log10.pval") {y = "log10.pval"}

  df2 = data.frame(x = df[,gsub("-", "", x)], y = df[,gsub("-", "", y)])
  p <- ggplot(df2)
  if(plottype == "x-y"){
    p <- p + geom_point(aes(x,y)) +theme_bw() +
      xlab(label = x)+
      ylab(label = y)
  }

  if(plottype == "histogram"){
    p <- p + geom_histogram(mapping = aes(x),bins = 50) +
      theme_bw() +
      xlab(label = x)
  }
  p
  return(p)

}
# statistics_plot(my_dep,"log10.padj","log10.padj",contrast = "KO_vs_WT",plottype = "x-y")



PTM_get_results <- function (dep) {
  library(tidyr)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  row_data$mean <- rowMeans(assay(dep), na.rm = TRUE)
  centered <- assay(dep) - row_data$mean
  centered <- data.frame(centered) %>% rownames_to_column() %>%
    gather(ID, val, -rowname) %>% left_join(., data.frame(colData(dep)),
                                            by = "ID")
  centered <- group_by(centered, rowname, condition) %>% summarize(val = mean(val,
                                                                              na.rm = TRUE)) %>% mutate(val = signif(val, digits = 3)) %>%
    spread(condition, val)
  colnames(centered)[2:ncol(centered)] <- paste(colnames(centered)[2:ncol(centered)],
                                                "_centered", sep = "")
  ratio <- as.data.frame(row_data) %>% column_to_rownames("name") %>%
    dplyr::select(ends_with("_diff")) %>% signif(., digits = 3) %>%
    rownames_to_column()
  colnames(ratio)[2:ncol(ratio)] <- gsub("_diff", "_ratio",
                                         colnames(ratio)[2:ncol(ratio)])
  df <- left_join(ratio, centered, by = "rowname")
  pval <- as.data.frame(row_data) %>% column_to_rownames("name") %>%
    dplyr::select(ends_with("p.val"), ends_with("p.adj"), ends_with("significant")) %>%
    rownames_to_column()
  pval[, grep("p.adj", colnames(pval))] <- pval[, grep("p.adj",
                                                       colnames(pval))] %>% signif(digits = 3)
  ids <- as.data.frame(row_data) %>% dplyr::select(name, ID)
  table <- left_join(ids, pval, by = c(name = "rowname"))
  table <- left_join(table, df, by = c(name = "rowname")) %>%
    arrange(desc(significant))
  protein.name = table$name %>% sapply(.,function(x){strsplit(x,"_")[[1]][1]})
  protein.ID = table$ID %>% sapply(.,function(x){strsplit(x,"_")[[1]][1]})
  sites = table$name %>% sapply(.,function(x){strsplit(x,"_")[[1]][2]})
  table2 = cbind(table[,1:2],protein.name,protein.ID,sites,table[,3:ncol(table)])
  rownames(table2) = rownames(table)
  return(table2)
}


PTM_get_table = function (results, type = c("centered", "contrast"))
{
  assertthat::assert_that(is.data.frame(results))
  type <- match.arg(type)
  if (length(grep("_ratio$", colnames(results))) < 1) {
    stop("'[contrast]_ratio' columns are not present in '",
         deparse(substitute(results)), "'.\nRun get_results() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_centered$", colnames(results))) < 1) {
    stop("'[contrast]_centered' columns are not present in '",
         deparse(substitute(results)), "'.\nRun get_results() to obtain the required columns",
         call. = FALSE)
  }
  significant <- results %>% filter(significant) %>% dplyr::select(-significant)
  if (type == "centered") {
    cols <- grep("_ratio", colnames(significant))
    table <- significant[, -cols]
    colnames(table)[c(1, 2)] <- c("PTM Name", "PTM ID")
    colnames(table)[grep("significant", colnames(table))] <- gsub("[.]",
                                                                  " - ", colnames(table)[grep("significant", colnames(table))])
    colnames(table) <- gsub("_centered", "", colnames(table)) %>%
      gsub("[_]", " ", .)
  }
  if (type == "contrast") {
    cols <- grep("_centered", colnames(significant))
    table <- significant[, -cols]
    colnames(table)[c(1, 2)] <- c("PTM Name", "PTM ID")
    colnames(table)[grep("significant", colnames(table))] <- gsub("[.]",
                                                                  " - ", colnames(table)[grep("significant", colnames(table))])
    colnames(table) <- gsub("_ratio", "", colnames(table)) %>%
      gsub("[_]", " ", .)
  }
  return(table)
}

get_dep_ptm_siglist = function (dep, alpha = 0.05, lfc = 1, threshold_method = c("intersect","curve"), curvature=0.6, x0_fold = 2)
{
  if (is.integer(alpha))
    alpha <- as.numeric(alpha)
  if (is.integer(lfc))
    lfc <- as.numeric(lfc)
  if (is.integer(curvature))
    curvature <- as.numeric(curvature)
  if (is.integer(x0_fold))
    x0_fold <- as.numeric(x0_fold)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                          length(lfc) == 1, is.numeric(curvature), length(curvature) == 1, is.numeric(x0_fold), length(x0_fold) == 1)
  # row_data <- rowData(dep, use.names = FALSE) %>% as.data.frame()
  res = get_results(dep)
  if (any(!c("name", "ID") %in% colnames(res))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(dep)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(res))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }

  if( length(threshold_method)!=1 | (!threshold_method %in% c("intersect","curve")) ){
    stop("threshold_method should be one of 'intersect'/'curve'",call. = F)
  }

  cols_p <- grep("_p.adj", colnames(res))
  cols_diff <- grep("_ratio", colnames(res))
  cols_sig <- grep("_significant", colnames(res))

  genelist_env = get_dep_ptm_genelist_Env()
  genelist_all_env = get_all_genelist_Env()
  rm(list = ls(envir = genelist_env),envir = genelist_env)
  # rm(list = ls(envir = genelist_all_env),envir = genelist_all_env)

  for(i in 1:length(cols_p)){
    tem = res[,cols_sig[i]]
    # tem = res[,cols_p[i]]<alpha & abs(res[,cols_diff[i]] )> lfc
    temp = res[which(tem),c(1:2,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp)>0){
      colnames(temp) = c("symbol","ID","L2FC","p.adj","significant")
      temp$significant = TRUE
      temp = temp[ , c(1, 3, 2, 4, 5)]
      if(threshold_method == "intersect") {
        assign(value = temp, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)
      } else if(threshold_method == "curve") {
        assign(value = temp, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Sig",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Sig",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)

      }
      #
      #       assign(value = temp, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
      #              envir = genelist_env)
      #       assign(value = temp, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
      #              envir = genelist_all_env)
    }

    tem_up = (res[,cols_sig[i]]) & (res[,cols_diff[i]] > 0)
    # tem_up = (res[,cols_p[i]]<alpha) & (res[,cols_diff[i]] > lfc)
    temp_up = res[which(tem_up),c(1:2,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp_up)>0){
      colnames(temp_up) = c("symbol","ID","L2FC","p.adj","significant")
      temp_up$significant = TRUE
      temp_up = temp_up[ , c(1, 3, 2, 4, 5)]

      if(threshold_method == "intersect") {
        assign(value = temp_up, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Up",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_up, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Up",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)

      } else if(threshold_method == "curve") {
        assign(value = temp_up, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Up",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_up, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Up",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)

      }

    }

    tem_dn = (res[,cols_sig[i]]) & (res[,cols_diff[i]] < 0)
    # tem_dn = (res[,cols_p[i]]<alpha) & (res[,cols_diff[i]] < -lfc)
    temp_dn = res[which(tem_dn),c(1:2,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp_dn)>0){
      colnames(temp_dn) = c("symbol","ID","L2FC","p.adj","significant")
      temp_dn$significant = TRUE
      temp_dn = temp_dn[ , c(1, 3, 2, 4, 5)]
      if(threshold_method == "intersect") {
        assign(value = temp_dn, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_dn, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)
      } else if(threshold_method == "curve") {
        assign(value = temp_dn, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_dn, x = paste("DEP_PTM_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)
      }

    }
  }
  return(ls(envir = get_dep_ptm_genelist_Env()))

}
