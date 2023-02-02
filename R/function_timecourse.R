
#\' @param algo Character(1), the clustering method, see \code{\link[TCseq]{timeclust}}.

#' Time-course clustering
#'
#' Expression pattern clustering for timecourse or multi-group experiment via
#' \code{\link[TCseq]{timeclust}()}
#'
#' @param x a SummarizedExperiment or DEGdata object.
#' @param ht_mat Data.frame or matrix, the expresiong assay. If x is provided, be reset to assay(x).
#' @param exp_design Data.frame, the experiment design. If x is provided, be reset to colData(x)
#' @param groupby Character(1), group by which column in experiment design.
#' @param group_order Character vector or factor, the order of groups in the \code{groupby} column.
#' @param k Integer(1), the number of clusters.
#' @param color Character(1), sets the color panel (from \pkg{RColorBrewer}).
#' It can be a palette name in "RdBu", "RdYlBu", "RdYlGn", "BrBG", "PiYG", "PRGn", "PuOr", "RdGy", "Spectral".
#' @param col_limit Numeric(1), sets the outer limits of the color scale.
#' @param row_font_size Integer(1), font size of row name
#' @param col_font_size Integer(1), font size of columns label.
#' @param heatmap_width Numeric(1)
#' @param heatmap_height Numeric(1)
#' @param seed Integer(1) or NULL, the random seed.
#'
#' @return A list contains:
#' - a cluster heatmap(Heatmap object from ComplexHeatmap),
#' - the cluster result for heatmap(data.frame),
#' - the output of \code{\link[TCseq]{timeclust}} and the seed.
#' @export
#'
#' @importFrom  TCseq timeclust clustCluster clustData clustMembership
#'
#' @examples
#' # Load sample
#' data(Silicosis_pg)
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols, mode = "delim", sep = "_")
#'
#' filt <- filter_se(se, thr = 0, fraction = 0.3, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.05)
#' diff <- test_diff(imputed,type = "control", control = "PBS")
#' dep <- add_rejections(diff)
#'
#' # Expression pattern cluster
#' TC <- get_tc_cluster(get_signicant(dep),
#'                      group_order = c("PBS","W2","W4","W6","W9","W10"), # set the group order
#'                      heatmap_width = 2.5, heatmap_height = 5
#' )
#' TC$ht # Heatmap
#' head(TC$res)
get_tc_cluster <- function(x, ht_mat, exp_design, groupby = "condition",
                           group_order = NULL,
                           # algo = "cm",
                           k = 6, dist = "euclidean",
                           color = c("RdBu", "RdYlBu", "RdYlGn", "BrBG", "PiYG", "PRGn", "PuOr", "RdGy", "Spectral"),
                           col_limit = 4, row_font_size = 5,
                           col_font_size = 5, heatmap_width = 3, heatmap_height = 5,
                           seed = NULL
                           # ,...
                           ){
  algo = "cm"
  color <- match.arg(color)
  cat("aa")
  if(hasArg(x) && !is.null(x)){
    ht_mat <- assay(x)
    exp_design <- colData(x)
  }
  cat("bb")

  print("start timeclust")
  if(is.null(group_order)) group_order = exp_design[,groupby] %>% unique()

  # ht_mat = ht_mat_save
  ht_mat = ht_mat[,which(exp_design[,groupby] %in% group_order)]
  ht_mat = ht_mat - rowMeans(ht_mat)

  exp_design2 = exp_design[which(exp_design[,groupby] %in% group_order),]
  ht_mat2 <- group_order %>% map_df(., function(x){
    rowMeans(ht_mat[,which(exp_design2[,groupby] == x)])
  }) %>%  t()

  colnames(ht_mat2) = group_order
  rownames(ht_mat2) = rownames(ht_mat)

  col_order <- group_order %>% sapply(., function(x){which(exp_design[,groupby] %in% x)}) %>% unlist
  col_order_save <<- col_order
  ht_mat = ht_mat[,col_order]

  if(is.null(seed)) seed <- sample(1:2^15, 1)
  set.seed(seed)

  cat(paste0("beging 'cm' cluster by dist = ", "'", dist, "'\n"))
  tca <- TCseq::timeclust(ht_mat2, algo = algo, k = k, standardize = TRUE, dist = dist)
  cluster1 <- TCseq::clustCluster(tca)
  tcadata <- TCseq::clustData(tca)

  Membershipdata <- TCseq::clustMembership(tca)
  mebership <- vector()
  for (i in 1:length(cluster1)) {
    mebership[i] <- Membershipdata[names(cluster1)[i],cluster1[i]]
  }
  names(mebership) <- names(cluster1)

  res_df <- cbind(name = rownames(ht_mat), Timecourse_cluster = cluster1, mebership = mebership, ht_mat) %>% as.data.frame()

  col_fun_line = circlize::colorRamp2(
    c(0.01, 0.5,0.9),
    c("#54c556","#d3ca48","#cd3634"))

  col_fun_heatmap = circlize::colorRamp2(
    seq(-col_limit, col_limit, (col_limit/5)),
    rev(RColorBrewer::brewer.pal(11, color)))

  ngroup = length(group_order)

  panel_fun2 = function(index, nm) {
    pushViewport(viewport(xscale = c(0,(2*ngroup+1)), yscale = c(-2,2)))
    grid.rect()

    gby = annotation_axis_grob(at = c(-1.7,0,1.7), labels = c(-1.7,0,1.7), labels_rot = 0,
                               side = "left", facing = "inside",gp = gpar(fontsize=5))
    # gby2 = annotation_axis_grob(at = 0,labels = "expression", labels_rot = 90,
    #                            side = "left", facing = "outside",gp = gpar(fontsize=5))
    grid.polyline(y=as.vector(t( (tcadata[index,]/4)+0.5 )),
                  x=c(rep(seq(from=2/(2*ngroup+2),to=(2*ngroup)/(2*ngroup+1),by=2/(2*ngroup+1) ), length(index))),
                  id=rep(1:length(index), times=rep(ncol(tcadata),times=length(index))),
                  gp=gpar(col=col_fun_line(mebership[index]), lwd=1,alpha = 0.5))

    if(cluster1[index[1]]==max(cluster1)){
      gbx = annotation_axis_grob(at = seq(from=2,to=(2*ngroup), by=2), labels = colnames(tcadata), labels_rot = 45,
                                 side = "bottom", facing = "outside",gp = gpar(fontsize= col_font_size ))
      grid.draw(gbx)
    }

    grid.draw(gby)
    # grid.draw(gby2)
    popViewport()
  }

  anno = anno_link(align_to = cluster1, which = "row", panel_fun = panel_fun2, side ="right",
                   size = unit((heatmap_height * 5/k) - 0.25, "cm"), gap = unit(0, "cm"), width = unit(heatmap_width*1.7, "cm"))

  ht <- ComplexHeatmap::Heatmap(ht_mat,
                                heatmap_width = unit(heatmap_width * 5, "cm"),heatmap_height = unit(heatmap_height * 5, "cm"),
                                # width = unit(5.5, "cm"),height = unit(20, "cm"),
                                col = col_fun_heatmap,
                                split = cluster1,
                                cluster_rows = F,
                                show_row_names = T,
                                row_names_side = "left",
                                row_dend_width = unit(10, "mm"),
                                cluster_row_slices = T,
                                cluster_column_slices =F,
                                row_gap = unit(1, "mm"),
                                #row_km = km,
                                cluster_columns =FALSE,
                                # column_labels =extractcolname(data),
                                heatmap_legend_param = list(
                                  title = "Centered value",
                                  title_position = "leftcenter-rot"
                                  # height = unit(500,"mm"),
                                ),
                                row_names_gp = gpar(fontsize = row_font_size),
                                column_names_gp = gpar(fontsize = col_font_size),
                                right_annotation = rowAnnotation(foo = anno)
  )

  print("finish heatmap")
  return(res_save <- list(res = res_df, ht = ht, mat = ht_mat, timeclust_res = tca, seed = seed))
}



