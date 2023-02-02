library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(limma)
library(colourpicker)
#require(tcltk)

coef_variation<-function(x){
  coef=sd(x)/mean(x)
}

#### Plot CVs

plot_cvs<-function(se) {

  ## backtransform data
  untransformed_intensity<- 2^(assay(se))
  exp_design<-colData(se)

### merge untransformed to exp design and calculate cvs

  cvs_group<- untransformed_intensity %>% data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("ID", "Intensity", -rowname) %>%
    dplyr::left_join(.,data.frame(exp_design), by="ID") %>%
    dplyr::group_by(rowname,condition) %>%
    dplyr::summarise(cvs=coef_variation(Intensity)) %>%
    dplyr::group_by(condition)%>%
    dplyr::mutate(condition_median=median(cvs))

p1 <-  ggplot(cvs_group, aes(cvs, color=condition, fill=condition)) +
    geom_histogram(alpha=.5, bins= 20, show.legend = FALSE) +
    facet_wrap(~condition) +
    geom_vline(aes(xintercept=condition_median, group=condition),color='grey40',
             linetype="dashed") +
    labs(title= 'Sample Coefficient of Variation', x="Coefficient of Variation", y="Count") +
    theme_DEP2() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold"))

p1 +geom_text(aes(x=max(cvs_group$cvs)-0.6,
                  y=max(ggplot_build(p1)$data[[1]]$ymax*1.1),
                  label=paste0("Median =",round(condition_median,2)*100,"%",by="")),
              show.legend = FALSE, size=4)

}




#manual: logical, if manual heatmap, default: FALSE
#manual_name: act when manual is TRUE, can be one or more contrast format "XXX_vs_XXX", default NULL. change to manual_contrast
#same_trend: same_trend:if the heatmap is the same trend, default TRUE; if the number of maual_name is 2(one from pro1, one from pro2),then the heatmap show: both significant and both trend(the two _vs_), if the number of maual_name is 4,6,8,10,...,then the heatmap show: same trend of significant(any significant), and same trend diff(log2(fc)) from the two part(_vs_), eg: two part is 1,1,1-2,2,2, same signicant can be T,T,T-T,T,T, T,F.F-T,F,F, F,T,F-F,T,F;same trend diff can be (---,---； -++，-++， "-" represent negative value, "+"represent positive and 0 value)
#if_mysplit: if order splices by yourself, default: FALSE
#mysplit: default: NULL, only act when if_mysplit is TRUE and kmeans is TRUE, it is numeric of 1:k, the numeric order is the splices order of heatmap(from top to bottom), eg: c(6,4,1,2,3,5)
#if_rowname_color: if set rowname color which Peptides == 1 to purple, Peptides == 2 to blue, default: FALSE
#if_chooseToshow: if the heatmap only to show proteins those you select only
#chooseToshow: act when if_chooseToshow is TRUE, the name of proteins those you want to show
# plot_heatmap <- function (dep, type = c("contrast", "centered"), manual = FALSE, manual_name = NULL, same_trend = FALSE,
#                           kmeans = FALSE, k = 6, color = "RdBu", col_limit = 6, indicate = NULL,
#                           clustering_distance = c("euclidean",
#                                                   "maximum", "manhattan", "canberra",
#                                                   "binary", "minkowski", "pearson", "spearman",
#                                                   "kendall", "gower"), row_font_size = 6, col_font_size = 10,
#                           plot = TRUE, if_mysplit = FALSE, mysplit = NULL, if_rowname_color = FALSE, if_chooseToshow = FALSE, chooseToshow = NULL, column_order = NULL, ...)
# {
#   library(ComplexHeatmap)
#   library(tidyr)
#   library(tibble)
#   if (is.integer(k))
#     k <- as.numeric(k)
#   if (is.integer(col_limit))
#     col_limit <- as.numeric(col_limit)
#   if (is.integer(row_font_size))
#     row_font_size <- as.numeric(row_font_size)
#   if (is.integer(col_font_size))
#     col_font_size <- as.numeric(col_font_size)
#   assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
#                           is.character(type), is.logical(kmeans), is.numeric(k),
#                           length(k) == 1, is.numeric(col_limit), length(col_limit) ==
#                             1, is.numeric(row_font_size), length(row_font_size) ==
#                             1, is.numeric(col_font_size), length(col_font_size) ==
#                             1, is.logical(plot), length(plot) == 1)
#
#   # set rowname of heatmap color, which Peptides == 1 to purple, ==2 to blue
#   if(if_rowname_color) {
#     rowData(dep)$row_name_color = "black"
#     rowData(dep)$row_name_color[which(rowData(dep)$Peptides == 1)] = "purple"
#     rowData(dep)$row_name_color[which(rowData(dep)$Peptides == 2)] = "blue"
#   }
#
#   ht_dep <<- dep
#   type <- match.arg(type)
#   clustering_distance <- match.arg(clustering_distance)
#   row_data <- rowData(dep, use.names = FALSE)
#   col_data <- colData(dep) %>% as.data.frame()
#
#
#
#   if (any(!c("label", "condition", "replicate") %in%
#           colnames(col_data))) {
#     stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
#                 deparse(substitute(dep)), "'"), call. = FALSE)
#   }
#   if (length(grep("_diff", colnames(row_data))) < 1) {
#     stop(paste0("'[contrast]_diff' columns are not present in '",
#                 deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."),
#          call. = FALSE)
#   }
#   if (!"significant" %in% colnames(row_data)) {
#     stop(paste0("'significant' column is not present in '",
#                 deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required column."),
#          call. = FALSE)
#   }
#   if (!is.null(indicate) & type == "contrast") {
#     warning("Heatmap annotation only applicable for type = 'centered'",
#             call. = FALSE)
#   }
#   if (!is.null(indicate) & type == "centered") {
#     ha1 <- get_annotation(dep, indicate)
#   } else {
#     ha1 <- NULL
#   }
#   #filtered <- dep[row_data$significant, ]
#   ## Filter for significant proteins only
#   library(stringr)
#   if(manual) {
#     #if it is both trend heatmap, a(data frame), b(check if both significant), you can add a if class(a),class(a) == "DFrame", set constrast muti select,if select two or more, heatmap is desplay both trend heatmap, you can changshi zhe zhong fang fa, tongli ratio, ind ye xugai,zai zai xiafang zhao xugai de difang, eg: all(c(1,0) < 0) | all(c(1,0)>=0)
#     #a = row_data[ , paste(manual_name, "_significant", sep = "")]
#     #b = apply(a, 1, all)
#     my_ind <- row_data[ , paste(manual_name, "_significant", sep = "")]
#     my = my_ind
#     if(class(my_ind) == "DFrame"){
#       if(same_trend){
#         my1 = apply(my_ind, 1, any)
#         #my_ind <- row_data[ , c(paste(manual_name, "_significant", sep = ""), "significant")]
#         my_ind = apply(my_ind, 1, function(i){
#           #all(i[1:((ncol(my_ind)-1)/2)] == i[((ncol(my_ind)-1 )/2 + 1):(ncol(my_ind)-1)]) & i[ncol(my_ind)]
#           all(i[1:((ncol(my_ind))/2)] == i[((ncol(my_ind) )/2 + 1):(ncol(my_ind))])
#         })
#         ratio <- row_data[ , paste(manual_name, "_diff", sep = "")]
#         ratio <- as.matrix(ratio)
#         ratio[which(ratio >= 0)] = 1
#         ratio[which(ratio < 0)] = -1
#
#         ratio <- apply(ratio, 1, function(myratio){
#           all(myratio[1:(ncol(ratio)/2)] == myratio[(ncol(ratio)/2 + 1):ncol(ratio)])
#         })
#
#         my_ind = intersect(which(my_ind),which(ratio))
#         my_ind  = intersect(my_ind, which(my1))
#
#       } else {
#         my_ind = apply(my_ind, 1, any)
#       }
#     }
#     filtered <- dep[my_ind , ]
#
#
#     #get index of manual_name
#     colname = colnames(assay(filtered))
#     # if(class(my) == "DFrame"){
#     #   ind = sort(unique(as.numeric(unlist(sapply(manual_name, simplify = FALSE, function(x){
#     #     which(str_detect(x, gsub("_[[:digit:]]*", "", colname)))
#     #   })))))
#     # } else {
#     #   # when this manual_name can only have one element
#     #   # #this is not common use: eg: when colnames(assay(dep)) = "Biotin_1" "Biotin_2" "Biotin_3" "H2A_1"    "H2A_2"    "H2A_3"    "uH2A_1"   "uH2A_2"   "uH2A_3" , set manual_name = "uH2A_vs_Biotin"
#     #   ind = which(str_detect(manual_name, gsub("_[[:digit:]]*", "", colname)))
#     #   #i = strsplit(manual_name, split = "_")[[1]][c(1,4)]# when DEP shiny app select parse from coloumn,this is c(1,4), else it is c(1,3), because parse from cols, treat_1 hui.bian.cheng treat__1
#     #   #ind = unlist(sapply(i, function(x){
#     #   #    which(str_detect(colname, x))
#     #   #    }))
#
#     # }
#     ind = which(gsub("_[[:digit:]]*", "", colname) %in% unlist(strsplit(manual_name, split = "_vs_")))
#   } else {
#     filtered <- dep[row_data$significant, ]
#   }
#
#   significant_save <<- row_data$significant
#
#   # Check for missing values
#   if(manual){
#     if(any(is.na(assay(filtered)[, ind]))) {
#       warning("Missing values in '", deparse(substitute(dep)), "'. ",
#               "Using clustering_distance = 'gower'",
#               call. = FALSE)
#       clustering_distance <- "gower"
#       obs_NA <- TRUE
#     } else {
#       obs_NA <- FALSE
#     }
#
#   } else{
#     if(any(is.na(assay(filtered)))) {
#       warning("Missing values in '", deparse(substitute(dep)), "'. ",
#               "Using clustering_distance = 'gower'",
#               call. = FALSE)
#       clustering_distance <- "gower"
#       obs_NA <- TRUE
#     } else {
#       obs_NA <- FALSE
#     }
#
#   }
#   # if(any(is.na(assay(filtered)))) {
#   #   warning("Missing values in '", deparse(substitute(dep)), "'. ",
#   #     "Using clustering_distance = 'gower'",
#   #     call. = FALSE)
#   #   clustering_distance <- "gower"
#   #   obs_NA <- TRUE
#   # } else {
#   #   obs_NA <- FALSE
#   # }
#
#   # Get centered intensity values ('centered')
#   if(type == "centered") {
#     if(manual){
#       rowData(filtered)$mean <- rowMeans(assay(filtered)[ , ind], na.rm = TRUE)
#       df <- assay(filtered)[ , ind] - rowData(filtered, use.names = FALSE)$mean
#     } else{
#       rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
#       df <- assay(filtered) - rowData(filtered, use.names = FALSE)$mean
#     }
#   }
#
#   filtered_save <<- filtered
#   # Get contrast fold changes ('contrast')
#   if(type == "contrast") {
#     df <- rowData(filtered, use.names = FALSE) %>%
#       data.frame() %>%
#       column_to_rownames(var = "name") %>%
#       dplyr::select(ends_with("_diff"))
#     colnames(df) <-
#       gsub("_diff", "", colnames(df)) %>%
#       gsub("_vs_", " vs ", .)
#     df <- as.matrix(df)
#     if(manual){
#       i = gsub("_vs_", " vs ", manual_name)
#       ii = as.data.frame(df[ , i])
#       colnames(ii) = i
#       df = as.matrix(ii)
#     }
#   }
#
#   # Facultative kmeans clustering
#   if(kmeans & obs_NA) {
#     warning("Cannot perform kmeans clustering with missing values",
#             call. = FALSE)
#     kmeans <- FALSE
#   }
#   if(kmeans & !obs_NA) {
#     set.seed(1)
#     df_kmeans <- kmeans(df, k)
#     if(type == "centered") {
#       # Order the k-means clusters according to the maximum fold change
#       # in all samples averaged over the proteins in the cluster
#       order <- data.frame(df) %>%
#         cbind(., cluster = df_kmeans$cluster) %>%
#         mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
#         group_by(cluster) %>%
#         summarize(index = sum(row)/n()) %>%
#         arrange(desc(index)) %>%
#         pull(cluster) %>%
#         match(seq_len(k), .)
#       df_kmeans$cluster <- order[df_kmeans$cluster]
#     }
#     if(type == "contrast") {
#       # Order the k-means clusters according to their average fold change
#       order <- data.frame(df) %>%
#         cbind(df, cluster = df_kmeans$cluster) %>%
#         gather(condition, diff, -cluster) %>%
#         group_by(cluster) %>%
#         summarize(row = mean(diff)) %>%
#         arrange(desc(row)) %>%
#         pull(cluster) %>%
#         match(seq_len(k), .)
#       df_kmeans$cluster <- order[df_kmeans$cluster]
#     }
#   }
#
#   # if(ncol(df) == 1) {
#   #   col_clust = FALSE
#   # } else {
#   #   col_clust = TRUE
#   # }
#   # if(nrow(df) == 1) {
#   #   row_clust = FALSE
#   # } else {
#   #   row_clust = TRUE
#   # }
#   if(clustering_distance == "gower") {
#     clustering_distance <- function(x) {
#       dist <- cluster::daisy(x, metric = "gower")
#       dist[is.na(dist)] <- max(dist, na.rm = TRUE)
#       return(dist)
#     }
#   }
#
#   # Legend info
#   legend <- ifelse(type == "contrast",
#                    "log2 Fold change",
#                    "log2 Centered intensity")
#
#   if(if_mysplit) {
#     cluster_row_slices = FALSE
#     cluster_column_slices = FALSE
#     split = if(kmeans) {factor(df_kmeans$cluster, levels = mysplit)} else {NULL}
#   } else {
#     cluster_row_slices = TRUE
#     cluster_column_slices = TRUE
#     split = if(kmeans) {df_kmeans$cluster} else {NULL}
#   }
#
#   #set rowname color
#   row_name_color = if(if_rowname_color) {rowData(filtered)$row_name_color} else {"black"}
#   if(if_chooseToshow) {
#     if(length(chooseToshow) >= 1) {
#       if(length(chooseToshow) == 1) {
#         df1 = df[match(chooseToshow, rownames(df)), ]
#         df1 = t(as.data.frame(df1))
#         rownames(df1) = rownames(df)[match(chooseToshow, rownames(df))]
#         df = df1
#       } else {
#         df = df[match(chooseToshow, rownames(df)), ]
#       }
#       row_name_color = if(if_rowname_color) {rowData(filtered)$row_name_color[match(chooseToshow, rownames(df))]} else {"black"}
#
#     } else {
#       return("at least one proteins should be selected")
#     }
#   }
#
#   ht_df <- df
#   if(!is.null(column_order)){
#     if(any(!column_order %in% colnames(df)))
#       return("at least one proteins should be selected")
#     df = df[,column_order]
#   }
#
#   # Heatmap
#   if(!if_chooseToshow) {
#     ht1 = ComplexHeatmap::Heatmap(df,
#                                   col = circlize::colorRamp2(
#                                     seq(-col_limit, col_limit, (col_limit/5)),
#                                     rev(RColorBrewer::brewer.pal(11, color))),
#                                   split = split,
#                                   cluster_row_slices = cluster_row_slices,
#                                   cluster_column_slices = cluster_column_slices,
#                                   # cluster_rows = col_clust,
#                                   # cluster_columns = row_clust,
#                                   row_names_side = "left",
#                                   column_names_side = "top",
#                                   clustering_distance_rows = clustering_distance,
#                                   clustering_distance_columns = clustering_distance,
#                                   heatmap_legend_param = list(color_bar = "continuous",
#                                                               legend_direction = "horizontal",
#                                                               legend_width = unit(5, "cm"),
#                                                               title_position = "lefttop"),
#                                   name = legend,
#                                   row_names_gp = gpar(fontsize = row_font_size, col = row_name_color),
#                                   column_names_gp = gpar(fontsize = col_font_size),
#                                   top_annotation = ha1,
#                                   ...)
#
#   } else {
#     ht1 = ComplexHeatmap::Heatmap(df,
#                                   col = circlize::colorRamp2(
#                                     seq(-col_limit, col_limit, (col_limit/5)),
#                                     rev(RColorBrewer::brewer.pal(11, color))),
#                                   # split = split,
#                                   # cluster_row_slices = cluster_row_slices,
#                                   # cluster_column_slices = cluster_column_slices,
#                                   # cluster_rows = col_clust,
#                                   # cluster_columns = row_clust,
#                                   row_names_side = "left",
#                                   column_names_side = "top",
#                                   clustering_distance_rows = clustering_distance,
#                                   clustering_distance_columns = clustering_distance,
#                                   heatmap_legend_param = list(color_bar = "continuous",
#                                                               legend_direction = "horizontal",
#                                                               legend_width = unit(5, "cm"),
#                                                               title_position = "lefttop"),
#                                   name = legend,
#                                   row_names_gp = gpar(fontsize = row_font_size, col = row_name_color),
#                                   column_names_gp = gpar(fontsize = col_font_size),
#                                   top_annotation = ha1,
#                                   ...)
#
#   }
#
#   if(plot) {
#     # Plot
#     draw(ht1, heatmap_legend_side = "top")
#   } else {
#     # Return data.frame
#     colnames(df) <- gsub(" ", "_", colnames(df))
#     df <- df[, unlist(column_order(ht1))]
#     if(kmeans) {
#       df <- cbind(df, k = df_kmeans$cluster)
#     }
#     return <- df[unlist(row_order(ht1)),]
#     data.frame(protein = row.names(return), return) %>%
#       mutate(order = row_number())
#   }
# }

# new add than DEP package
# add_values: if add values in the heatmap
# value_size: the size of added values, act when add_values is TRUE
# digits: integer indicating the number of decimal places
plot_dist <- function (dep, significant = TRUE, pal = "YlOrRd", pal_rev = TRUE,
                       indicate = NULL, font_size = 12, plot = TRUE, add_values = FALSE, value_size = 10, digits = 2, ...)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(significant), length(significant) == 1, is.character(pal),
                          length(pal) == 1, is.logical(pal_rev), length(pal_rev) ==
                            1, is.numeric(font_size), length(font_size) == 1,
                          is.logical(plot), length(plot) == 1)
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_dist() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- colData(dep) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", deparse(substitute(dep)),
           ".\nValid columns are: '", paste(columns,
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- colData(dep) %>% data.frame() %>% select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1)
        cols <- c("black")
      if (length(var) == 2)
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Pastel1")
      if (length(var) >= 7)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }
  else {
    ha1 <- NULL
  }
  if (significant) {
    if (!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
      stop("'significant' column is not present in '",
           deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required column",
           call. = FALSE)
    }
    dep <- dep[rowData(dep, use.names = FALSE)$significant,
               ]
  }
  dist_mat <- cluster::daisy(t(assay(dep)), metric = "gower") %>%
    as.matrix()
  max <- max(dist_mat)
  ht1 = Heatmap(dist_mat, col = circlize::colorRamp2(seq(0,
                                                         max, ((max)/7)), if (pal_rev) {
                                                           rev(RColorBrewer::brewer.pal(8, pal))
                                                         }
                                                     else {
                                                       RColorBrewer::brewer.pal(8, pal)
                                                     }), heatmap_legend_param = list(color_bar = "continuous",
                                                                                     legend_direction = "horizontal", legend_width = unit(5,
                                                                                                                                          "cm"), title_position = "topcenter"),
                name = "Gower's distance", column_names_gp = gpar(fontsize = font_size),
                row_names_gp = gpar(fontsize = font_size), top_annotation = ha1, cell_fun = if(add_values) {function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf(paste("%.", digits, "f", sep = ""), dist_mat[i, j]), x, y, gp = gpar(fontsize = value_size))}} else {NULL},
                ...)
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(dist_mat)
    return(df)
  }
}

# new add than DEP package
# add_values: if add values in the heatmap
# value_size: the size of added values, act when add_values is TRUE
# digits: integer indicating the number of decimal places
plot_cor <- function (dep, significant = TRUE, lower = -1, upper = 1, pal = "PRGn",
                      pal_rev = FALSE, indicate = NULL, font_size = 12, plot = TRUE,  add_values = FALSE, value_size = 10, digits = 2,
                      ...)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(significant), length(significant) == 1, is.numeric(lower),
                          length(lower) == 1, is.numeric(upper), length(upper) ==
                            1, is.character(pal), length(pal) == 1, is.logical(pal_rev),
                          length(pal_rev) == 1, is.numeric(font_size), length(font_size) ==
                            1, is.logical(plot), length(plot) == 1)
  if (!(lower >= -1 & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid\n         Run plot_cor() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_cor() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  if (any(is.na(assay(dep)))) {
    stop("Missing values in '", deparse(substitute(dep)),
         "'. Use plot_dist() instead")
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- colData(dep) %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", deparse(substitute(dep)),
           ".\nValid columns are: '", paste(columns,
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- colData(dep) %>% data.frame() %>% select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1)
        cols <- c("black")
      if (length(var) == 2)
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Pastel1")
      if (length(var) >= 7)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }
  else {
    ha1 <- NULL
  }
  if (significant) {
    if (!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
      stop("'significant' column is not present in '",
           deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required column",
           call. = FALSE)
    }
    dep <- dep[rowData(dep, use.names = FALSE)$significant,
               ]
  }
  cor_mat <- cor(assay(dep))
  ht1 = Heatmap(cor_mat, col = circlize::colorRamp2(seq(lower,
                                                        upper, ((upper - lower)/7)), if (pal_rev) {
                                                          rev(RColorBrewer::brewer.pal(8, pal))
                                                        }
                                                    else {
                                                      RColorBrewer::brewer.pal(8, pal)
                                                    }), heatmap_legend_param = list(color_bar = "continuous",
                                                                                    legend_direction = "horizontal", legend_width = unit(5,
                                                                                                                                         "cm"), title_position = "topcenter"),
                name = "Pearson correlation", column_names_gp = gpar(fontsize = font_size),
                row_names_gp = gpar(fontsize = font_size), top_annotation = ha1, cell_fun = if(add_values) {function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf(paste("%.", digits, "f", sep = ""), cor_mat[i, j]), x, y, gp = gpar(fontsize = value_size))}} else {NULL},
                ...)
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(cor_mat)
    return(df)
  }
}



plot_missval <- function (se, ...)
{
    assertthat::assert_that(inherits(se, "SummarizedExperiment"))
    se_assay <- assay(se)
    if (!any(is.na(se_assay))) {
        stop("No missing values in '", deparse(substitute(se)),
            "'", call. = FALSE)
    }
    df <- se_assay %>% data.frame(.)
    missval <- df[apply(df, 1, function(x) any(is.na(x))), ]
    missval <- ifelse(is.na(missval), 0, 1)
    ht2 = ComplexHeatmap::Heatmap(missval, col = c("white", "black"),
        column_names_side = "top", show_row_names = FALSE,
        show_column_names = TRUE, name = "Missing values pattern",
        column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(at = c(0,
            1), labels = c("Missing value", "Valid value")), ...)
    draw(ht2, heatmap_legend_side = "top")
}
############  new add parameters
#same_width: if axis same width from 0, default TRUE
#my_breaks: if set x breaks by yourself, default FALSE
#mybreaks: only my_breaks = TRUE act, the numeric vector, eg: seq(-10,10,by = 5), c(-2, -1, 0, 1, 2)
plot_volcano <- function (dep, contrast, label_size = 3, add_names = TRUE, adjusted = FALSE,
    plot = TRUE, same_width = TRUE, my_breaks = FALSE, mybreaks = NULL)
{
    if (is.integer(label_size))
        label_size <- as.numeric(label_size)
    assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
        is.character(contrast), length(contrast) == 1, is.numeric(label_size),
        length(label_size) == 1, is.logical(add_names), length(add_names) ==
            1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot),
        length(plot) == 1)
    row_data <- rowData(dep, use.names = FALSE)
    if (any(!c("name", "ID") %in% colnames(row_data))) {
        stop(paste0("'name' and/or 'ID' columns are not present in '",
            deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."),
            call. = FALSE)
    }
    if (length(grep("_p.adj|_diff", colnames(row_data))) <
        1) {
        stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
            deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."),
            call. = FALSE)
    }
    if (length(grep("_significant", colnames(row_data))) <
        1) {
        stop(paste0("'[contrast]_significant' columns are not present in '",
            deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."),
            call. = FALSE)
    }
    if (length(grep(paste(contrast, "_diff", sep = ""),
        colnames(row_data))) == 0) {
        valid_cntrsts <- row_data %>% data.frame() %>% dplyr::select(ends_with("_diff")) %>%
            colnames(.) %>% gsub("_diff", "", .)
        valid_cntrsts_msg <- paste0("Valid contrasts are: '",
            paste0(valid_cntrsts, collapse = "', '"), "'")
        stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
            valid_cntrsts_msg, call. = FALSE)
    }
    #when name have Disadvantages eg: when contrast = c("H2A_vs_Biotin","uH2A_vs_Biotin"), when set contrast = c("H2A_vs_Biotin"), can find two cols,and then error, the same problem with p_values and signif
    #diff <- grep(paste(contrast, "_diff", sep = ""),
    #    colnames(row_data))
    diff <- match(paste(contrast, "_diff", sep = ""),
        colnames(row_data))
    if (adjusted) {
        p_values <- match(paste(contrast, "_p.adj", sep = ""),
            colnames(row_data))
    }
    else {
        p_values <- match(paste(contrast, "_p.val", sep = ""),
            colnames(row_data))
    }
    signif <- match(paste(contrast, "_significant", sep = ""),
        colnames(row_data))
    df <- data.frame(x = row_data[, diff], y = -log10(row_data[,
        p_values]), significant = row_data[, signif], name = row_data$name) %>%
        dplyr::filter(!is.na(significant)) %>% arrange(significant)
    name1 <- gsub("_vs_.*", "", contrast)
    name2 <- gsub(".*_vs_", "", contrast)
    xlimit <- ceiling(max(c(abs(min(df$x)), abs(max(df$x)))))
    p <- ggplot(df, aes(x, y)) + geom_vline(xintercept = 0) +
        geom_point(aes(col = significant)) +
        geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1,
            0), vjust = c(-1, -1), label = c(name1, name2), size = 5, fontface = "bold")) +
        #scale_x_continuous(limits = c(-xlimit, xlimit), breaks = seq(-15, 15, by = 5)) +
        labs(title = contrast,
        x = expression(log[2] ~ "Fold change")) + theme_DEP1() +
        theme(legend.position = "none") + scale_color_manual(values = c(`TRUE` = "black",
        `FALSE` = "grey"))

    if (add_names) {
        p <- p + ggrepel::geom_text_repel(data = dplyr::filter(df, significant),
            aes(label = name), size = label_size, box.padding = unit(0.1,
                "lines"), point.padding = unit(0.1, "lines"),
            segment.size = 0.5)

    }
    if (adjusted) {
        p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
    }
    else {
        p <- p + labs(y = expression(-log[10] ~ "P-value"))

    }
    if(same_width){
        p <- p + scale_x_continuous(limits = c(-xlimit, xlimit))
    }
    if(same_width & my_breaks){
        p <- p + scale_x_continuous(limits = c(-xlimit, xlimit), breaks = mybreaks)
    }
    if (plot) {
        return(p)
    }
    else {
        df <- df %>% dplyr::select(name, x, y, significant) %>% arrange(desc(x))
        colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change",
            "p_value_-log10")
        if (adjusted) {
            colnames(df)[3] <- "adjusted_p_value_-log10"
        }
        return(df)
    }
}




## for RNAseq differential analysis
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(BiocParallel)
# library(apeglm) ## shrink
library(DESeq2)
# library(org.Hs.eg.db)
# library(AnnotationDbi)
library(dplyr)
## library(IHW)
# library(ashr)
# /**
#  * [pca plot]
#  * @param  {[type]} data       [row gene , column sample]
#  * @param  {[type]} coldata    [dataframe contain label, ID,condition and replicate,eg:]
#  # DataFrame with 12 rows and 4 columns
# #              label          ID   condition replicate
# #        <character> <character> <character> <numeric>
# # Ubi4_1      Ubi4_1      Ubi4_1        Ubi4         1
# # Ubi4_2      Ubi4_2      Ubi4_2        Ubi4         2
# # Ubi4_3      Ubi4_3      Ubi4_3        Ubi4         3
# # Ubi6_1      Ubi6_1      Ubi6_1        Ubi6         1
# # Ubi6_2      Ubi6_2      Ubi6_2        Ubi6         2
# # ...            ...         ...         ...       ...
# # Ctrl_2      Ctrl_2      Ctrl_2        Ctrl         2
# # Ctrl_3      Ctrl_3      Ctrl_3        Ctrl         3
# # Ubi1_1      Ubi1_1      Ubi1_1        Ubi1         1
# # Ubi1_2      Ubi1_2      Ubi1_2        Ubi1         2
# # Ubi1_3      Ubi1_3      Ubi1_3        Ubi1         3
#  * @param  {Number} x          [Integer(1), Sets the principle component to plot on the x-axis]
#  * @param  {Number} y          [Integer(1), Sets the principle component to plot on the y-axis]
#  * @param  {[type]} indicate   [Character, Sets the color, shape and facet_wrap of the plot based on columns from the experimental coldata]
#  * @param  {[type]} label      [Logical, Whether or not to add sample labels]
#  * @param  {Number} n          [Integer(1), Sets the number of top variable proteins to consider]
#  * @param  {Number} point_size [Integer(1), Sets the size of the points]
#  * @param  {Number} label_size [Integer(1), Sets the size of the labels]
#  * @param  {[type]} plot       [Logical(1), If TRUE (default) the PCA plot is produced. Otherwise (if FALSE), the data which the PCA plot is based on are returned]
#  * @param  {[type]} if_square  [if xlim == ylim, default FALSE]
#  * @return {[type]}            [description]
#  */
plot_my_pca <- function (data, coldata, x = 1, y = 2, indicate = c("condition", "replicate"), label = FALSE, n = 500, point_size = 4,
                         label_size = 3, plot = TRUE, if_square = FALSE)
{
  library(ggplot2)
  library(DEP)
  library(dplyr)
  if (is.integer(x))
    x <- as.numeric(x)
  if (is.integer(y))
    y <- as.numeric(y)
  if (is.integer(n))
    n <- as.numeric(n)
  if (is.integer(point_size))
    point_size <- as.numeric(point_size)
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)

  assertthat::assert_that(
    is.numeric(x), length(x) == 1, is.numeric(y), length(y) ==
      1, is.numeric(n), length(n) == 1, is.character(indicate),
    is.logical(label), is.numeric(point_size), length(point_size) ==
      1, is.numeric(label_size), length(label_size) ==
      1, is.logical(plot), length(plot) == 1)

  if (x > ncol(data) | y > ncol(data)) {
    stop(paste0("'x' and/or 'y' arguments are not valid\n",
                "Run plot_my_pca() with 'x' and 'y' <= ", ncol(data),
                "."), call. = FALSE)
  }
  if (n > nrow(data)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_my_pca() with 'n' <= ",
                nrow(data), "."), call. = FALSE)
  }
  # for shiny indicate can select the third indicate
  mycolData <- coldata
  mycolData$Condition = mycolData$condition
  mycolData$Replicate = mycolData$replicate

  columns <- colnames(mycolData)
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_my_pca() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ", deparse(substitute(data)),
                  ".\nValid columns are: '", paste(columns,
                                                   collapse = "', '"), "'."), call. = FALSE)
    }
  }
  var <- apply(data, 1, sd)
  df <- data[order(var, decreasing = TRUE)[seq_len(n)],
             ]
  pca <- prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame() %>% rownames_to_column() %>%
    left_join(., mycolData, by = c(rowname = "ID"))
  percent <- round(100 * pca$sdev^2/sum(pca$sdev^2), 1)
  for (feat in indicate) {
    pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }
  limit = unlist(pca_df[,c(paste0("PC", x), paste0("PC", y))])
  p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC", y)))) +
  labs(title = paste0("PCA plot - top ", n, " variable genes"), x = paste0("PC", x, ": ", percent[x], "%"), y = paste0("PC", y, ": ", percent[y], "%")) +
    coord_fixed() +
    theme_DEP1()
  if(if_square){
    p <- p +
      scale_x_continuous(limits = range(limit)) +
      scale_y_continuous(limits = range(limit))
  }

  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]]),
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]],
                            shape = pca_df[[indicate[2]]]), size = point_size) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = pca_df[[indicate[1]]],
                            shape = pca_df[[indicate[2]]]), size = point_size) +
      facet_wrap(~pca_df[[indicate[3]]]) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    return(p)
  }
  else {
    df <- pca_df %>% dplyr::select(rowname, paste0("PC", c(x,
                                                           y)), match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

# /**
#  * [generates a distance matrix heatmap using the Gowers distance]
#  * @param  {[type]}    data      [row gene , column sample]
#  * @param  {[type]}    coldata   [dataframe contain label, ID,condition and replicate]
#  * @param  {String}    pal       [Character(1), Sets the color panel (from RColorBrewer)]
#  * @param  {[type]}    pal_rev   [Logical(1), Whether or not to invert the color palette]
#  * @param  {[type]}    indicate  [Character, Sets additional annotation on the top of the heatmap based on columns from the experimental design (coldata)]
#  * @param  {Number}    font_size [Integer(1), Sets the size of the labels]
#  * @param  {[type]}    plot      [Logical(1), If TRUE (default) the distance matrix plot is produced. Otherwise (if FALSE), the data which the distance matrix plot is based on are returned]
#  * @param  {...[type]}           [Additional arguments for Heatmap function as depicted in Heatmap]
#  * @return {[type]}              [description]
#  */
#  new add
# add_values: if add values in the heatmap
# value_size: the size of added values, act when add_values is TRUE
# # digits: integer indicating the number of decimal places
plot_my_dist <- function (data, coldata, pal = "YlOrRd", pal_rev = TRUE,
                          indicate = NULL, font_size = 12, plot = TRUE, add_values = FALSE, value_size = 10, digits = 2, ...)
{
  library(ComplexHeatmap)
  library(RColorBrewer)
  assertthat::assert_that(
    is.character(pal),
    length(pal) == 1, is.logical(pal_rev), length(pal_rev) ==
      1, is.numeric(font_size), length(font_size) == 1,
    is.logical(plot), length(plot) == 1)
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_my_dist() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- coldata
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", deparse(substitute(dep)),
           ".\nValid columns are: '", paste(columns,
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- coldata %>% dplyr::select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1)
        cols <- c("black")
      if (length(var) == 2)
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Pastel1")
      if (length(var) >= 7)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }
  else {
    ha1 <- NULL
  }
  # if (significant) {
  #     if (!"significant" %in% colnames(coldata)) {
  #         stop("'significant' column is not present in '",
  #             deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required column",
  #             call. = FALSE)
  #     }
  #     data <- data[data$significant,
  #         ]
  # }
  dist_mat <- cluster::daisy(t(data), metric = "gower") %>%
    as.matrix()
  max <- max(dist_mat)
  ht1 = Heatmap(dist_mat, col = circlize::colorRamp2(seq(0, max, ((max)/7)), if (pal_rev) {rev(RColorBrewer::brewer.pal(8, pal))}
                                                     else {RColorBrewer::brewer.pal(8, pal)
                                                     }), heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal", legend_width = unit(5, "cm"), title_position = "topcenter"), name = "Gower's distance", column_names_gp = gpar(fontsize = font_size),
                row_names_gp = gpar(fontsize = font_size), top_annotation = ha1, cell_fun = if(add_values) {function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf(paste("%.", digits, "f", sep = ""), dist_mat[i, j]), x, y, gp = gpar(fontsize = value_size))}} else {NULL},
                ...)
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(dist_mat)
    return(df)
  }
}

# /**
#  * [generates a Pearson correlation matrix]
#  * @param  {[type]}    data      [row gene , column sample]
#  * @param  {[type]}    coldata   [dataframe contain label, ID,condition and replicate]
#  * @param  {Number}    lower     [Integer(1), Sets the lower limit of the color scale]
#  * @param  {Number}    upper     [Integer(1), Sets the upper limit of the color scale]
#  * @param  {String}    pal       [Character(1), Sets the color panel (from RColorBrewer)]
#  * @param  {[type]}    pal_rev   [Logical(1), Whether or not to invert the color palette]
#  * @param  {[type]}    indicate  [Character, Sets additional annotation on the top of the heatmap based on columns from the experimental design (coldata)]
#  * @param  {Number}    font_size [Integer(1), Sets the size of the labels]
#  * @param  {[type]}    plot      [Logical(1), If TRUE (default) the correlation matrix plot is produced. Otherwise (if FALSE), the data which the correlation matrix plot is based on are returned]
#  * @param  {...[type]}           [Additional arguments for Heatmap function as depicted in Heatmap]
#  * @return {[type]}              [description]
#  */
#  #  new add
# add_values: if add values in the heatmap
# value_size: the size of added values, act when add_values is TRUE
# digits: integer indicating the number of decimal places

plot_my_cor <- function (data, coldata, lower = -1, upper = 1, pal = "PRGn",
                         pal_rev = FALSE, indicate = NULL, font_size = 12, plot = TRUE, add_values = FALSE, value_size = 10, digits = 2,
                         ...)
{
  library(ComplexHeatmap)
  library(RColorBrewer)
  assertthat::assert_that(
    is.numeric(lower),
    length(lower) == 1, is.numeric(upper), length(upper) ==
      1, is.character(pal), length(pal) == 1, is.logical(pal_rev),
    length(pal_rev) == 1, is.numeric(font_size), length(font_size) ==
      1, is.logical(plot), length(plot) == 1)
  if (!(lower >= -1 & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid\n         Run plot_my_cor() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_my_cor() with one of the following 'pal' options: ",
         paste(pals$rowname, collapse = "', '"), "'",
         call. = FALSE)
  }
  if (any(is.na(data))) {
    stop("Missing values in '", deparse(substitute(data)),
         "'. Use plot_dist() instead")
  }
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- coldata %>% as.data.frame()
    columns <- colnames(col_data)
    if (any(!indicate %in% columns)) {
      stop("'", paste0(indicate, collapse = "' and/or '"),
           "' column(s) is/are not present in ", deparse(substitute(dep)),
           ".\nValid columns are: '", paste(columns,
                                            collapse = "', '"), "'.", call. = FALSE)
    }
    anno <- coldata %>% data.frame() %>% dplyr::select(indicate)
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var = anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1)
        cols <- c("black")
      if (length(var) == 2)
        cols <- c("orangered", "cornflowerblue")
      if (length(var) < 7 & length(var) > 2)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Pastel1")
      if (length(var) >= 7)
        cols <- RColorBrewer::brewer.pal(length(var),
                                         "Set3")
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 = HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }
  else {
    ha1 <- NULL
  }
  # if (significant) {
  #     if (!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
  #         stop("'significant' column is not present in '",
  #             deparse(substitute(dep)), "'\nRun add_rejections() to obtain the required column",
  #             call. = FALSE)
  #     }
  #     data <- dep[data$significant,
  #         ]
  # }
  cor_mat <- cor(data)
  ht1 = Heatmap(cor_mat, col = circlize::colorRamp2(seq(lower,
                                                        upper, ((upper - lower)/7)), if (pal_rev) {
                                                          rev(RColorBrewer::brewer.pal(8, pal))
                                                        }
                                                    else {
                                                      RColorBrewer::brewer.pal(8, pal)
                                                    }), heatmap_legend_param = list(color_bar = "continuous",
                                                                                    legend_direction = "horizontal", legend_width = unit(5,
                                                                                                                                         "cm"), title_position = "topcenter"),
                name = "Pearson correlation", column_names_gp = gpar(fontsize = font_size),
                row_names_gp = gpar(fontsize = font_size), top_annotation = ha1, cell_fun = if(add_values) {function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf(paste("%.", digits, "f", sep = ""), cor_mat[i, j]), x, y, gp = gpar(fontsize = value_size))}} else {NULL},
                ...)
  if (plot) {
    draw(ht1, heatmap_legend_side = "top")
  }
  else {
    df <- as.data.frame(cor_mat)
    return(df)
  }
}


########## 如果symbol列在res_obj则intgenes用symbol列，若不在，则用rownames(res_obj)列
#' @param res_obj A \code{\link{DESeqResults}} object,the result of function get_res or the result of adding SYMBOL following by function get_res
#' @param contrast your contrast,formed as A_vs_B
#' @param FDR Numeric value, the significance level for thresholding adjusted p-values
#' @param log2fc  Numeric value, the significance level for thresholding log2 fold change
#' @param point_alpha Alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param sig_color Color to use to mark differentially expressed genes. Defaults to red
#' e.g. HGNC-based gene symbols. Optional
#' @param hlines The y coordinate (in absolute value) where to draw horizontal lines,
#' optional
#' @param title A title for the plot, optional
#' @param xlab X axis label, defaults to "mean of normalized counts - log10 scale"
#' @param ylim Vector of two numeric values, Y axis limits to restrict the view
#' @param add_rug Logical, whether to add rug plots in the margins
#' @param intgenes Vector of genes of interest. Gene symbols if a \code{symbol}
#' column is provided in \code{res_obj}, or else the identifiers specified in the
#' row names
#' @param intgenes_color The color to use to mark the genes on the main plot.
#' @param labels_intgenes Logical, whether to add the gene identifiers/names close
#' to the marked plots
#' @param labels_repel Logical, whether to use \code{geom_text_repel} for placing the
#' labels on the features to mark
#'
#' @return An object created by \code{ggplot}
#' @export
plot_ma <- function(res_obj,
                    contrast,
                    FDR = 0.05,
                    log2fc = 2,
                    point_alpha = 0.2,
                    sig_color = "red3",
                    hlines = NULL,
                    title = contrast,
                    xlab = "mean of normalized counts - log10 scale",
                    ylim = NULL,
                    add_rug = FALSE,
                    intgenes = NULL,
                    intgenes_color = "steelblue4",
                    labels_intgenes = TRUE,
                    labels_repel = TRUE) {
  library(ggrepel)
  ma_df <- data.frame(
    mean = res_obj[ , paste(contrast, "baseMean", sep = "_")],
    lfc = res_obj[ , paste(contrast, "log2FoldChange", sep = "_")],
    padj = res_obj[ , paste(contrast, "padj", sep = "_")],
    isDE = ifelse(is.na(res_obj[ , paste(contrast, "padj", sep = "_")]), FALSE, res_obj[ , paste(contrast, "padj", sep = "_")] < FDR & abs(res_obj[ , paste(contrast, "log2FoldChange", sep = "_")]) >= log2fc),
    ID = rownames(res_obj)
  )

  ma_df <- ma_df[ma_df$mean > 0, ]

  ma_df$logmean <- log10(ma_df$mean) # TO ALLOW FOR BRUSHING!!


  p <- ggplot(ma_df, aes_string(x = "logmean", y = "lfc", colour = "isDE"))

  if (!is.null(hlines)) {
    p <- p + geom_hline(aes(yintercept = hlines), col = "steelblue4", alpha = 1) +
      geom_hline(aes(yintercept = -hlines), col = "steelblue4", alpha = 1)
  }
  p <- p + geom_hline(aes(yintercept = 0), col = "red3", alpha = 1)

  p <- p + xlab(xlab) + ylab("log2 fold change")

  p <- p + geom_point(alpha = point_alpha)
  p <- p + scale_colour_manual(
    name = paste0("FDR = ", FDR, "\n", "log2 fc = ", log2fc),
    values = c("TRUE" = sig_color, "FALSE" = "gray21"),
    labels = c("TRUE" = "DE", "FALSE" = "nonDE")
  )

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }

  if (!is.null(intgenes)) {

    # now here for the symbol
    res_df <- as.data.frame(res_obj)
    res_df$logmean <- log10(res_df[, paste(contrast, "baseMean", sep = "_")])

    if ("symbol" %in% colnames(res_df)) {
      # use the gene names
      df_intgenes <- res_df[res_df$symbol %in% intgenes, ]
      df_intgenes$myids <- df_intgenes$symbol
    } else {
      # use whatever is there as id
      df_intgenes <- res_df[rownames(res_df) %in% intgenes, ]
      df_intgenes$myids <- rownames(df_intgenes)
    }

    # df_intgenes <- res_df[res_df$symbol %in% intgenes,]
    p <- p + geom_point(data = df_intgenes, aes_string("logmean", paste(contrast, "log2FoldChange", sep = "_")), color = intgenes_color, size = 4)

    if (labels_intgenes) {
      if (labels_repel) {
        p <- p + geom_text_repel(
          data = df_intgenes, aes_string("logmean", paste(contrast, "log2FoldChange", sep = "_"), label = "myids"),
          color = intgenes_color, size = 5
        )
      } else {
        p <- p + geom_text(
          data = df_intgenes, aes_string("logmean", paste(contrast, "log2FoldChange", sep = "_"), label = "myids"),
          color = intgenes_color, size = 5, hjust = 0.25, vjust = -0.75
        )
      }
    }
  }

  if (add_rug) {
    p <- p + geom_rug(alpha = 0.3)
  }

  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  p
}


# /**
#  * [generates a heatmap for RNAseq]
#  * @param  {[type]}    data                [data matrix of ntd,rownames(data) == rownames(res),res is the result of function get_res,default: rownames(data) is row label]
#  * @param  {[type]}    row_data            [the result of function get_res or the result of adding SYMBOL following by function get_res ]
#  * @param  {[type]}    type                [one of c("raw", "centered")]
#  * @param  {[type]}    manual              [logical, if manual heatmap, default: FALSE]
#  * @param  {[type]}    manual_name         [act when manual is TRUE, can be one or more contrast format "XXX_vs_XXX", default NULL]
#  * @param  {[type]}    kmeans              [Logical(1), Whether or not to perform k-means clustering]
#  * @param  {Number}    k                   [Integer(1), Sets the number of k-means clusters]
#  * @param  {Number}    col_limit           [Integer(1), Sets the outer limits of the color scale]
#  * @param  {[type]}    clustering_distance ["euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall" or "gower", Function used to calculate clustering distance (for proteins and samples). Based on Heatmap and daisy]
#  * @param  {Number}    row_font_size       [Integer(1), Sets the size of row labels]
#  * @param  {Number}    col_font_size       [Integer(1), Sets the size of column labels]
#  * @param  {[type]}    plot                [Logical(1), If TRUE (default) the heatmap is produced. Otherwise (if FALSE), the data which the heatmap is based on are returned]
#  * @param  {[type]}    if_mysplit          [if order splices by yourself, default: FALSE]
#  * @param  {[type]}    mysplit             [default: NULL, only act when if_mysplit is TRUE and kmeans is TRUE, it is numeric of 1:k, the numeric order is the splices order of heatmap(from top to bottom), eg: c(6,4,1,2,3,5)]
#  * @param  {...[type]}                     [Additional arguments for Heatmap function as depicted in Heatmap]
#  * @return {[type]}                        [description]
#  */
#if_chooseToshow: if the heatmap only to show proteins those you select only
#chooseToshow: act when if_chooseToshow is TRUE, the name of proteins those you want to show

plot_heatmap_rnaseq <- function (data, row_data, type = c("raw", "centered"), manual = FALSE, manual_name = NULL, kmeans = FALSE, k = 6, color = "RdBu", col_limit = 6, clustering_distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall", "gower"), row_font_size = 6, col_font_size = 10,plot = TRUE, if_mysplit = FALSE, mysplit = NULL, if_chooseToshow = FALSE, chooseToshow = NULL, ...)
{
  library(ComplexHeatmap)
  library(tidyr)
  library(tibble)
  if (is.integer(k))
    k <- as.numeric(k)
  if (is.integer(col_limit))
    col_limit <- as.numeric(col_limit)
  if (is.integer(row_font_size))
    row_font_size <- as.numeric(row_font_size)
  if (is.integer(col_font_size))
    col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(is.character(type), is.logical(kmeans), is.numeric(k),
                          length(k) == 1, is.numeric(col_limit), length(col_limit) ==
                            1, is.numeric(row_font_size), length(row_font_size) ==
                            1, is.numeric(col_font_size), length(col_font_size) ==
                            1, is.logical(plot), length(plot) == 1, all(rownames(data) == rownames(row_data)) | all(na.omit(rownames(data)) == na.omit(row_data$symbol)))

  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  row_data <- row_data

  library(stringr)
  if(manual) {
    my_ind <- row_data[ , paste(manual_name, "_significant", sep = "")]
    my = my_ind
    if(class(my_ind) == "data.frame"){
        my_ind = apply(my_ind, 1, any)

    }
    filtered <- data[my_ind , ]

    #get index of manual_name
    colname = colnames(filtered)

    ind = which(gsub("_[[:digit:]]*", "", colname) %in% unlist(strsplit(manual_name, split = "_vs_")))
  } else {
    filtered <- data[row_data$significant, ]
  }

  data = filtered

  # Check for missing values
  if(manual){
    if(any(is.na(filtered[, ind]))) {
      warning("Missing values in '", deparse(substitute(filtered)), "'. ",
              "Using clustering_distance = 'gower'",
              call. = FALSE)
      clustering_distance <- "gower"
      obs_NA <- TRUE
    } else {
      obs_NA <- FALSE
    }

  } else{
    if(any(is.na(filtered))) {
      warning("Missing values in '", deparse(substitute(filtered)), "'. ",
              "Using clustering_distance = 'gower'",
              call. = FALSE)
      clustering_distance <- "gower"
      obs_NA <- TRUE
    } else {
      obs_NA <- FALSE
    }

  }

  # Get centered intensity values ('centered')
  if(type == "centered") {
    if(manual){
      filtered = filtered[ , ind]
      df <- filtered - rowMeans(filtered)
    } else {
      df <- filtered - rowMeans(filtered)
    }
  } else {
    if(manual) {
      df = filtered[ , ind]
    } else {
      df <- filtered
    }
  }

  # Facultative kmeans clustering
  if(kmeans & obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
            call. = FALSE)
    kmeans <- FALSE
  }
  if(kmeans & !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    # Order the k-means clusters according to the maximum fold change
    # in all samples averaged over the proteins in the cluster
    order <- data.frame(df) %>%
      cbind(., cluster = df_kmeans$cluster) %>%
      mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
      group_by(cluster) %>%
      summarize(index = sum(row)/n()) %>%
      arrange(desc(index)) %>%
      pull(cluster) %>%
      match(seq_len(k), .)
    df_kmeans$cluster <- order[df_kmeans$cluster]
  }


  if(clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }

  # Legend info
  legend <- ifelse(type == "raw",
                   "normalized count",
                   "centered normalized count")

  if(if_mysplit) {
    cluster_row_slices = FALSE
    cluster_column_slices = FALSE
    split = if(kmeans) {factor(df_kmeans$cluster, levels = mysplit)} else {NULL}
  } else {
    cluster_row_slices = TRUE
    cluster_column_slices = TRUE
    split = if(kmeans) {df_kmeans$cluster} else {NULL}
  }

  if(type == "raw"){
    col = colorRampPalette(rev(brewer.pal(n = 9, name =color)))(100)
  } else {
    col = circlize::colorRamp2(
      seq(-col_limit, col_limit, (col_limit/5)),
      rev(RColorBrewer::brewer.pal(11, color)))
  }

  if(if_chooseToshow) {
    if(length(chooseToshow) >= 1) {
      if(length(chooseToshow) == 1) {
        df1 = df[match(chooseToshow, rownames(df)), ]
        df1 = t(as.data.frame(df1))
        rownames(df1) = rownames(df)[match(chooseToshow, rownames(df))]
        df = df1
      } else {
        df = df[match(chooseToshow, rownames(df)), ]
      }

    } else {
      return("at least one proteins should be selected")
    }
  }

  # Heatmap
  if(!if_chooseToshow) {
    ht1 = ComplexHeatmap::Heatmap(df,
                                  col = col,
                                  split = split,
                                  cluster_row_slices = cluster_row_slices,
                                  cluster_column_slices = cluster_column_slices,
                                  # cluster_rows = col_clust,
                                  # cluster_columns = row_clust,
                                  row_names_side = "left",
                                  column_names_side = "top",
                                  clustering_distance_rows = clustering_distance,
                                  clustering_distance_columns = clustering_distance,
                                  heatmap_legend_param = list(color_bar = "continuous",
                                                              legend_direction = "horizontal",
                                                              legend_width = unit(5, "cm"),
                                                              title_position = "lefttop"),
                                  name = legend,
                                  row_names_gp = gpar(fontsize = row_font_size),
                                  column_names_gp = gpar(fontsize = col_font_size),
                                  ...)

  } else {
    ht1 = ComplexHeatmap::Heatmap(df,
                                  col = col,
                                  # split = split,
                                  # cluster_row_slices = cluster_row_slices,
                                  # cluster_column_slices = cluster_column_slices,
                                  # cluster_rows = col_clust,
                                  # cluster_columns = row_clust,
                                  row_names_side = "left",
                                  column_names_side = "top",
                                  clustering_distance_rows = clustering_distance,
                                  clustering_distance_columns = clustering_distance,
                                  heatmap_legend_param = list(color_bar = "continuous",
                                                              legend_direction = "horizontal",
                                                              legend_width = unit(5, "cm"),
                                                              title_position = "lefttop"),
                                  name = legend,
                                  row_names_gp = gpar(fontsize = row_font_size),
                                  column_names_gp = gpar(fontsize = col_font_size),
                                  ...)
  }

  if(plot) {
    # Plot
    draw(ht1, heatmap_legend_side = "top")
  } else {
    # Return data.frame
    colnames(df) <- gsub(" ", "_", colnames(df))
    df <- df[, unlist(column_order(ht1))]
    if(kmeans) {
      df <- cbind(df, k = df_kmeans$cluster)
    }
    return <- df[unlist(row_order(ht1)),]
    data.frame(protein = row.names(return), return) %>%
      mutate(order = row_number())
  }
}

# /**
#  * [generates a volcano plot for a specified contrast, if column symbol is in res, the label is symbol, else is rownames(res)]
#  * @param  {[type]} res             [the result of diffentical analysis, no any row ordered]
#  * @param  {[type]} contrast        [Character(1), Specifies the contrast to plot, form as "trt_vs_untrt"]
#  * @param  {[type]} adjusted        [Logical(1), Whether or not to use adjusted p values]
#  * @param  {String} labelWay        [the way of adding label, one of c("all significant", "up","down","selected proteins"), default: "all significant"]
#  * @param  {Number} showNum         [the number of the added labels corresponding to the labelWay, default: 20, when showNum = 0, there is no label added, and when it is more than the labelWay, it is the nrow(labelWay), eg: you can set it 10000 in order to show all label of the labelWay]
#  * @param  {[type]} chooseTolabel   [the protein names of your choose , when labelWay == "selected proteins" act , default: NULL]
#  * @param  {Number} fontSize        [the font size of the label, default: 7]
#  * @param  {Number} dotsize         [the dot size of the label, default: 0.8]
#  * @param  {[type]} same_width      [if axis same width from 0, default TRUE]
#  * @param  {[type]} plot            [If TRUE (default) the volcano plot is produced. Otherwise (if FALSE), the data which the volcano plot is based on are returned.]
#  * @param  {[type]} fcCutoff        [log2 fold change cutoff]
#  * @param  {Number} linetype        [linetype for hline or vline, default:4]
#  * @param  {Number} adjpCutoff      [adj pvalue cutoff]
#  * @param  {[type]} label.rectangle [TRUE or FALSE, whether rectangle around label]
#  * @param  {Number} stroke          [label point lwd, default : 0.6]
#  * @return {[type]}                 [description]
#  */
plot_volcano_rnaseq <- function(res, contrast, adjusted = FALSE ,labelWay = "up",showNum = 20, chooseTolabel = NULL,fontSize = 9, dotsize = 2, same_width = TRUE, plot = TRUE, fcCutoff = log2(2), linetype = 4, adjpCutoff = 0.05, label.rectangle = FALSE, stroke = 0.6, down_color = "#0099CC", stable_color = "#BEBEBE", up_color = "#CD8500") {
  library(ggpubr)
  library(ggthemes)
  library(DEP)

  row_data <- res

  diff <- match(paste(contrast, "log2FoldChange", sep = "_"), colnames(row_data))
  p_values <- match(paste(contrast, "pvalue", sep = "_"), colnames(row_data))
  fdr <- match(paste(contrast, "padj", sep = "_"), colnames(row_data))

  if("symbol" %in% colnames(row_data)) {
    name = row_data$symbol
  } else {
    name = rownames(row_data)
  }


  if(adjusted) {
    df <- data.frame(x = row_data[, diff], y = -log10(row_data[, fdr]), fdr = row_data[, fdr], name = name)
  } else {
    df <- data.frame(x = row_data[, diff], y = -log10(row_data[, p_values]), fdr = row_data[, fdr], name = name)
  }
  df$significant = df$fdr < adjpCutoff & abs(df$x) >= fcCutoff
  df <- df %>% filter(!is.na(significant)) %>% arrange(significant)

  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  xlimit <- ceiling(max(c(abs(min(df$x)), abs(max(df$x)))))


  test <- df
  test$change = ifelse((test$significant),
                       ifelse(test$x > 0 ,'up','down'),
                       'stable')

  if(labelWay == "selected proteins"){
    labels <- chooseTolabel
    for_label <- test[which(test$name %in% labels),]
    if(is.null(labels) | nrow(for_label) == 0) return(print("please choose at least one protein"))
    ########################################################到这里
  }else if(labelWay == "all significant"){
    labelProtein <- test[which(test$change!='stable'),]
    labelProtein <- labelProtein[order(labelProtein$y,decreasing = T),]
    labels <- labelProtein$name[c(0:min(nrow(labelProtein),showNum))]
    labels <- na.omit(labels)
    for_label <- test[which(test$name %in% labels),]
  }else{
    labelProtein <-test[which(test$change == labelWay),]
    labelProtein <- labelProtein[order(labelProtein$y,decreasing = T),]
    labels <- labelProtein$name[c(0:min(nrow(labelProtein),showNum))]
    labels <- na.omit(labels)
    for_label <- test[which(test$name %in% labels),]
  }

  if(adjusted) {
    ylab = expression(-log[10] ~ "Adjusted p-value")
  } else {
    ylab = expression(-log[10] ~ "P-value")
  }
  volcanoPlot <-  ggscatter(test,
                            x = "x",
                            y = "y",
                            title = contrast,
                            xlab = expression(log[2] ~ "Fold change"),
                            ylab = ylab,
                            shape = 16,#20
                            color = "change",
                            size = dotsize,
                            label = "name",
                            font.label = c(fontSize, "plain") ,
                            repel = T,
                            palette = c("down" = down_color, 'stable' = stable_color, "up" = up_color) ,# #CC9933,#FC4E07,lightcoral,#0099CC,#00AFBB
                            label.select = labels,
                            label.rectangle = label.rectangle) +
    theme_base() +
    #adjust pvalue volcano in DEP is very ugly, Basically not used, so no need to add the limit line based on padj
    # geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col="grey",lwd=0.6) +
    geom_vline(xintercept=c(-fcCutoff,fcCutoff),linetype=linetype, col=c(down_color, up_color),lwd=0.6) +
    theme(plot.title=element_text(hjust=0.5)) +
    theme_DEP1() +
    theme(legend.position = "none") +
    # geom_vline(xintercept = 0) +
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1,
                                                                                      0), vjust = c(-1, -1), label = c(name1, name2), size = 5, fontface = "bold"))+
    # scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expand_scale(mult = 0, add = c(0, 0))) +
    geom_point(aes(fill = change), size = if(labelWay == "selected proteins"){dotsize + 1} else {dotsize}, shape = 21, color = "black", stroke = stroke, position = "identity", data = for_label)

  if(same_width){
    volcanoPlot <- volcanoPlot + scale_x_continuous(limits = c(-xlimit, xlimit))
  }

  if(adjusted){#for adjust pvalue
    volcanoPlot <- volcanoPlot + geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col=stable_color,lwd=0.6)
  } else {#for pvalue
          volcanoPlot <- volcanoPlot + geom_hline(yintercept=min(test$y[!test$change == "stable"]),linetype=linetype, col=stable_color,lwd=0.6)
      #     volcanoPlot <- volcanoPlot + geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col="grey",lwd=0.6)
      }
  volcanoPlot
  # else {#for pvalue
  #     # volcanoPlot <- volcanoPlot + geom_hline(yintercept=min(test$y[!test$change == "stable"]),linetype=linetype, col="grey",lwd=0.6)
  #     volcanoPlot <- volcanoPlot + geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col="grey",lwd=0.6)
  # }

  if (plot) {
    return(volcanoPlot)
  }
  else {
    test <- test %>% dplyr::select(name, x, y, significant, change) %>% arrange(desc(x))
    colnames(test)[c(1, 2, 3)] <- c("gene", "log2_fold_change",
                                    "p_value_-log10")
    if (adjusted) {
      colnames(test)[3] <- "adjusted_p_value_-log10"
    }
    return(test)
  }
}

# /**
#  * [get differential analysis result of a specified contrast]
#  * @param  {[type]} dds                  [the dds object returned from function DEseq]
#  * @param  {[type]} coldata                  [the exdesign]
#  * @param  {[type]} type             ["control", "all" or "manual", The type of contrasts that will be tested. This can be all possible pairwise comparisons ("all"), limited to the comparisons versus the control ("control"), or manually defined contrasts ("manual").]
#  * @param  {[type]} control                  [The condition to which contrasts are generated if type = "control" (a control condition would be most appropriate)]
#  * @param  {[type]} test                  [The contrasts that will be tested if type = "manual". These should be formatted as "SampleA_vs_SampleB" or c("SampleA_vs_SampleC", "SampleB_vs_SampleC")]
#  * @param  {[type]} contrast_upon                  [the contrast to bulid upon]
#  * @param  {[type]} filter_ihw           [whether to use ihw to perform independent filtering and p-value adjustment, default FALSE]
#  * @param  {[type]} independentFiltering [logical, whether independent filtering should be applied automatically, default TRUE]
#  * @param  {[type]} lfcshark             [whether to use lfcshark, default FALSE]
#  * @return {[type]}                      [description]
#  */
get_res <- function(dds, coldata, type = c("control", "all", "manual"), control = NULL, test = NULL, contrast_upon = "condition", filter_ihw = FALSE, independentFiltering = TRUE, lfcshark = FALSE) {
  type <- match.arg(type)
  conditions <- as.character(unique(coldata[ , contrast_upon]))

  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                    collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""),
                   cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control, "- ", sep = " "), "", .) %>% paste(" - ", control, sep = "")
      }
    }
  }
  if (type == "control") {
    if (is.null(control))
      stop("Please select control")
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = " - ")
  }
  if (type == "manual") {
    if (is.null(test)) {
      stop("Please select one or more contrasts that you want to test")
    }
    assertthat::assert_that(is.character(test))
    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("Valid contrasts should be like this: conditionA_vs_conditionB")
    }
    cntrst <- gsub("_vs_", " - ", test)
  }
  message("extracted results: ", paste(gsub(" - ",
                                            "_vs_", cntrst), collapse = ", "))
  getres <- function(dds, contrast, filter_ihw = FALSE, independentFiltering = TRUE, lfcshark = FALSE) {
    if(!filter_ihw) {
      res = results(dds, contrast = contrast, independentFiltering = independentFiltering)
    } else {
      res = results(dds, contrast = contrast, independentFiltering = independentFiltering, filterFun = ihw)
    }
    if(lfcshark){
      res = lfcShrink(dds = dds, contrast = contrast, res = res, type = "ashr")
    }
    return(res)
  }

  res_df <- sapply(cntrst, function(i){
    re <- getres(dds = dds, contrast = c(contrast_upon, gsub(" - .*", "", i), gsub(".* - ", "", i)), filter_ihw = filter_ihw, independentFiltering = independentFiltering, lfcshark = lfcshark)

  })
  names(res_df) = gsub(" - ","_vs_",names(res_df))
  res_df = do.call(cbind, res_df)
  names(res_df) = gsub("\\.", "_", names(res_df))

  test_res <- new("DEGdata",dds)
  test_res@test_result = res_df
  return(test_res)
}

# /**
#  * [Mark significant proteins]
#  * @param  {[type]} diff  [the result of function get_res or the result of adding SYMBOL following by function get_res]
#  * @param  {Number} alpha Sets the threshold for the adjusted P value.]
#  * @param  {Number} lfc   [Sets the threshold for the log2 fold change]
#  * @return {[type]}       [description]
#  */
add_rejection <- function (diff, alpha = 0.05, lfc = 1) {
  if (is.integer(alpha))
    alpha <- as.numeric(alpha)
  if (is.integer(lfc))
    lfc <- as.numeric(lfc)
  assertthat::assert_that(is.numeric(alpha), length(alpha) == 1,
                          is.numeric(lfc), length(lfc) == 1,
                          class(diff) =="DEGdata")
  # row_data <- as.data.frame(diff)
  diff_table <- diff@test_result
  row_data <- as.data.frame(diff_table)

  if (length(grep("_padj|_log2FoldChange", colnames(row_data))) <
      1) {
    stop("'[contrast]_log2FoldChange' and/or '[contrast]_padj' columns are not present in '",
         deparse(substitute(diff_table)), "'\nRun get_res() to obtain the required columns",
         call. = FALSE)
  }
  cols_p <- grep("_padj", colnames(row_data))
  cols_diff <- grep("_log2FoldChange", colnames(row_data))
  if (length(cols_p) == 1) {
    diff_table$significant <- row_data[, cols_p] <= alpha &
      abs(row_data[, cols_diff]) >= lfc
    diff_table$significant[which(is.na(diff_table$significant))] = FALSE
    diff_table$contrast_significant <- diff_table$significant
    colnames(diff_table)[ncol(diff_table)] <- gsub("padj","significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df, significant = apply(sign_df,
                                                  1, function(x) any(x)))
    colnames(sign_df) <- gsub("_padj", "_significant", colnames(sign_df))

    diff_table <- merge(diff_table, sign_df, by = "row.names")
    rownames(diff_table) = diff_table$Row.names
    diff_table = diff_table[,-1]
  }
  # diff$significant[which(is.na(diff$significant))] = FALSE
  diff@test_result <- diff_table
  return(diff)
}

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

#### function for protein LFQ imputation data and RNA normalized count data heatmap together
# /**
#  * [description]
#  * @param  {[type]} rna               [the normalized count matrix for RNAseq, the first column is gene, and other columns are your samples end with "_" + digit, eg: name,  A_1, A_2, A_3, B_1, B_2, B_3]
#  * @param  {[type]} pro               [the protein matrix eg: the log2, normalized and imputation protein matrix from R package DEP]
#  * @param  {[type]} gene              [data frame, the gene list that you want to plot heatmap, the column name is gene]
#  * @param  {[type]} Type              [Choose the genes to show. It can be one of 'all', 'both', 'at least one'. 'all' represents genes that all you uploaded; 'both' represents genes that both have RNAseq data and protein data; 'at least one' represents genes that at least have one of RNAseq and protein data]
#  * @param  {Number} row_font_size     [the row font size]
#  * @param  {Number} col_font_size     [the column font size]
#  * @param  {Number} column_title_size [the title of the heatmap]
#  * @param  {Number} Numcol_pro        [the number color limit of protein data]
#  * @param  {Number} Numcol_rna        [the number color limit of rnaseq data]
#  * @param  {[type]} cluster_columns   [if cluster columns]
#  * @return {[type]}                   [description]
#  */
#if_chooseToshow: if the heatmap only to show proteins those you select only
#chooseToshow: act when if_chooseToshow is TRUE, the name of proteins those you want to show
plot_pro_rna_heatmap <- function(rna, pro, gene, Type = "all", row_font_size = 6, col_font_size = 8, column_title_size = 10, Numcol_pro = 2, Numcol_rna = 2, cluster_columns = F, color = "RdBu", if_chooseToshow = FALSE, chooseToshow = NULL) {
  library(dplyr)
  library(circlize)
  library(RColorBrewer)
  library(ComplexHeatmap)
  ##### in order to know which gene no rna or pro data, and set gene color to grey which have no data
  pro$name_pro = pro$name
  rna$name_rna = rna$name

  if(Type == "all genes") {
    gene = gene
  }

  if(Type == "both") {
    gene = data.frame(name = intersect(rna$name[na.omit(match(gene$name, rna$name))], pro$name[na.omit(match(gene$name, pro$name))]), stringsAsFactors = F)
  }

  if(Type == "at least one") {
    gene = data.frame(name = union(rna$name[na.omit(match(gene$name, rna$name))], pro$name[na.omit(match(gene$name, pro$name))]), stringsAsFactors = F)
  }

  df_rna = rna[match(gene$name, rna$name),]
  df_rna$gene = gene$name
  df <- merge(df_rna, pro, by.x = "gene", by.y = "name", all.x = TRUE)
  #rna, pro set rowname color, no data to grey
  df$row_color_RNA = "black"
  df$row_color_RNA[is.na(df$name_rna)] = "grey"

  df$row_color_pro = "black"
  df$row_color_pro[is.na(df$name_pro)] = "grey"
  colnames(df) = gsub("\\.y$", "", gsub("\\.x$", "", colnames(df)))

  dat = df[ , c(match("gene", colnames(df)), match("name_rna", colnames(df)), match("row_color_RNA", colnames(df)), 3:(match("name_rna", colnames(df)) -1), match("name_pro", colnames(df)), match("row_color_pro", colnames(df)), (match("name_rna", colnames(df)) +1):(match("name_pro", colnames(df)) -1))]
  names(dat) = gsub("\\.1","",colnames(dat))


  if(if_chooseToshow) {
    if(length(chooseToshow) >= 1) {
      if(length(chooseToshow) == 1) {
        df1 = dat[match(chooseToshow, dat$gene), ]
        df1 = as.data.frame(df1)
        rownames(df1) = dat$gene[match(chooseToshow, dat$gene)]
        dat = df1
      } else {
        dat = dat[match(chooseToshow, dat$gene), ]
      }

    } else {
      return("at least one proteins should be selected")
    }
  }

  RNA = dat[ , match("row_color_RNA", colnames(dat)) : (match("name_pro", colnames(dat)) - 1)]
  PRO = dat[ , match("row_color_pro", colnames(dat)) : ncol(dat)]
  rownames(RNA) = dat$gene
  rownames(PRO) = dat$gene


  ################### PRO heatmap
  dat_pro_NA = PRO[apply(PRO[ , -1], 1, function(i){all(is.na(i))}), ]
  dat_pro = PRO[apply(PRO[ , -1], 1, function(i){!all(is.na(i))}), ]
  dat_pro_heat = dat_pro[ , -1] - rowMeans(dat_pro[ , -1])

  ht1=ComplexHeatmap::Heatmap(dat_pro_heat, col = colorRamp2(seq(-Numcol_pro, Numcol_pro, (Numcol_pro/5)), rev(brewer.pal(11, color))),
              na_col = "grey96",row_names_gp = gpar(fontsize = row_font_size),
              column_names_gp = gpar(fontsize = col_font_size),
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T,
              heatmap_legend_param = list(color_bar = "continuous",title_position = "topleft",legend_width = unit(5, "cm")),
              name = "Protein", cluster_columns = cluster_columns, row_km = 1, column_title = "Protein centered Enrichment (Log2)", column_title_gp = gpar(fontsize = column_title_size))

  rw = row_order(ht1)
  dat_pro = dat_pro[rw, ]
  dat_pro = rbind(dat_pro, dat_pro_NA)
  row_name_color_pro = dat_pro$row_color_pro
  dat_pro = dat_pro[ ,-1]
  dat_pro = dat_pro - rowMeans(dat_pro)

  ht1=ComplexHeatmap::Heatmap(dat_pro, col = colorRamp2(seq(-Numcol_pro, Numcol_pro, (Numcol_pro/5)), rev(brewer.pal(11, color))),
              na_col = "grey96",row_names_gp = gpar(fontsize = row_font_size, col = row_name_color_pro),
              column_names_gp = gpar(fontsize = col_font_size),
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T,
              heatmap_legend_param = list(color_bar = "continuous",title_position = "topleft",legend_width = unit(5, "cm")),
              name = "Protein", cluster_columns = cluster_columns, row_km = 1, column_title = "Protein centered Enrichment (Log2)", column_title_gp = gpar(fontsize = column_title_size), cluster_rows = FALSE)


  ht_pro <- draw(ht1, heatmap_legend_side = c("right"))

  ################### RNA heatmap
  dat_rna = RNA[match(rownames(dat_pro), rownames(RNA)), ]
  row_name_color_rna = dat_rna$row_color_RNA
  dat_rna = dat_rna[ ,-1]
  dat_rna = dat_rna - rowMeans(dat_rna)

  ht2=ComplexHeatmap::Heatmap(dat_rna, col = colorRamp2(seq(-Numcol_rna, Numcol_rna, (Numcol_rna/5)), rev(brewer.pal(11, color))),
              na_col = "grey96",row_names_gp = gpar(fontsize = row_font_size, col = row_name_color_rna),
              column_names_gp = gpar(fontsize = col_font_size),
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T,
              heatmap_legend_param = list(color_bar = "continuous",title_position = "topleft"), name = "RNA", cluster_columns = cluster_columns, row_km = 1,cluster_rows = F, column_title = "RNA centered normalized count", column_title_gp = gpar(fontsize = column_title_size))#width = unit(6, "cm")

  ht_rna <- draw(ht2, heatmap_legend_side = c("right"))

  ht_lis = ht1 + ht2
  hhh = draw(ht_lis, ht_gap = unit(0.4, "cm"), main_heatmap = "Protein", auto_adjust = FALSE,legend_title_position = "top",heatmap_legend_side = "right")

  return(lis = list(ht_pro = ht_pro, ht_rna = ht_rna, ht_merge = hhh, dat = dat))
}

#### for PR-Heatmap Connected to the whole app
# /**
#  * [description]
#  * @param  {[type]} rna [the full_results for downloading from DEG-RNAseq panel]
#  * @return {[type]}     [description]
#  */
prepare_rna <- function(rna) {
  names(rna)[1] = "name"
  names(rna) = paste("rna", names(rna), sep = "_")
  return(rna)
}

# /**
#  * [description]
#  * @param  {[type]} pro [the full_dataset for downloading from DEP-LFQ panel]
#  * @return {[type]}     [description]
#  */
prepare_pro <- function(pro) {
  index_1 = which(gsub("_[[:digit:]]*", "", names(pro)) %in% unique(unlist(strsplit(gsub("_significant", "", grep("_significant", names(pro), value = T)), split = "_vs_"))))
  index_2 = grep("_vs_|significant", names(pro))
  pro = pro[ , c(1, index_1, index_2)]
  names(pro) = paste("pro", names(pro), sep = "_")
  return(pro)
}

#### function for protein LFQ imputation data and RNA normalized count data heatmap together
# /**
#  * [description]
#  * @param  {[type]} rna               [the result from function prepare_rna]
#  * @param  {[type]} pro               [the result from function prepare_pro]
#  * @param  {[type]} manual            [if use manual contrast]
#  * @param  {[type]} contrast          [the contrast, act when manual is TRUE]
#  * @param  {[type]} Type              [Choose the genes to show. based on union RNAseq significant and protein significant, it can be one of 'both' 'both significant', and 'at least one' . both represents genes that both have RNAseq data and protein data; at least one represents genes that at least have one of RNAseq and protein data, both significant represents both RNAseq data and protein data are significant]
#  * @param  {Number} row_font_size     [the row font size]
#  * @param  {Number} col_font_size     [the column font size]
#  * @param  {Number} column_title_size [the title of the heatmap]
#  * @param  {Number} Numcol_pro        [the number color limit of protein data]
#  * @param  {Number} Numcol_rna        [the number color limit of rnaseq data]
#  * @param  {[type]} cluster_columns   [if cluster columns]
#  * @return {[type]}                   [description]
#  */
#if_chooseToshow: if the heatmap only to show proteins those you select only
#chooseToshow: act when if_chooseToshow is TRUE, the name of proteins those you want to show
plot_pro_rna_heatmap_for_whole <- function(rna, pro, manual = FALSE, contrast = NULL, Type = "at least one", row_font_size = 6, col_font_size = 8, column_title_size = 10, Numcol_pro = 2, Numcol_rna = 2, cluster_columns = F, color = "RdBu", if_chooseToshow = FALSE, chooseToshow = NULL) {
  library(dplyr)
  library(circlize)
  library(RColorBrewer)
  library(ComplexHeatmap)

  if(manual) {
  gene_rna = rna %>% dplyr::filter(get(paste("rna_", contrast, "_significant", sep = ""))) %>% dplyr::select(rna_name)
  gene_pro = pro %>% dplyr::filter(get(paste("pro_", contrast, "_significant", sep = ""))) %>% dplyr::select(pro_name)
  if(Type == "at least one"){
    gene = data.frame(name = na.omit(union(gene_rna$rna_name, gene_pro$pro_name)))
  } else {
    if(Type == "both significant") {
      gene = data.frame(name = na.omit(intersect(gene_rna$rna_name, gene_pro$pro_name)))
      } else {
        if(Type == "both") {
          gene = na.omit(union(gene_rna$rna_name, gene_pro$pro_name))
          gene = data.frame(name = intersect(rna$rna_name[na.omit(match(gene, rna$rna_name))], pro$pro_name[na.omit(match(gene, pro$pro_name))]), stringsAsFactors = F)
        }
      }
  }
}


  ###### any significant
  if(!manual) {
    gene_rna = rna %>% dplyr::filter(rna_significant) %>% dplyr::select(rna_name)
    gene_pro = pro %>% dplyr::filter(pro_significant) %>% dplyr::select(pro_name)
    if(Type == "at least one") {
      gene = data.frame(name = na.omit(union(gene_rna$rna_name, gene_pro$pro_name)))
    } else {
      if(Type == "both significant") {
        gene = data.frame(name = na.omit(intersect(gene_rna$rna_name, gene_pro$pro_name)))
        } else {
          if(Type == "both") {
            gene = na.omit(union(gene_rna$rna_name, gene_pro$pro_name))
            gene = data.frame(name = intersect(rna$rna_name[na.omit(match(gene, rna$rna_name))], pro$pro_name[na.omit(match(gene, pro$pro_name))]), stringsAsFactors = F)
          }
        }
    }
  }

  ##### in order to know which gene no rna or pro data, and set gene color to grey which have no data
  pro$name_pro = pro$pro_name
  rna$name_rna = rna$rna_name

  df_rna = rna[match(gene$name, rna$rna_name),]
  df_rna$gene = gene$name
  df <- merge(df_rna, pro, by.x = "gene", by.y = "pro_name", all.x = TRUE)

  ###### for contrast
  #rna, pro set rowname color, no data to grey
  if(manual) {
    df$row_color_RNA = "black"
    df$row_color_RNA[df[ , paste("rna_", contrast, "_significant", sep = "")]] = "red"
    df$row_color_RNA[is.na(df$name_rna)] = "grey"

    df$row_color_pro = "black"
    df$row_color_pro[df[ , paste("pro_", contrast, "_significant", sep = "")]] = "red"
    df$row_color_pro[is.na(df$name_pro)] = "grey"

    colnames(df) = gsub("rna_", "", gsub("pro_", "", colnames(df)))

    df = df[ , -grep("_vs_|significant", names(df))]
    names(df) = gsub("\\.1","",colnames(df))
    df = df[ , which(gsub("_\\d*$", "", colnames(df)) %in% c("gene", "name_rna", "name_pro", "row_color_RNA", "row_color_pro", unlist(strsplit(contrast, split = "_vs_"))))]
    names(df) = gsub("\\.1","",colnames(df))
    dat = df[ , c(match("gene", colnames(df)), match("name_rna", colnames(df)), match("row_color_RNA", colnames(df)), 2:(match("name_rna", colnames(df)) -1), match("name_pro", colnames(df)), match("row_color_pro", colnames(df)), (match("name_rna", colnames(df)) +1):(match("name_pro", colnames(df)) -1))]
    names(dat) = gsub("\\.1","",colnames(dat))
  }



  ##### for any significant
  #rna, pro set rowname color, no data to grey
  if(!manual) {
    df$row_color_RNA = "black"
    df$row_color_RNA[df[ , "rna_significant"]] = "red"
    df$row_color_RNA[is.na(df$name_rna)] = "grey"

    df$row_color_pro = "black"
    df$row_color_pro[df[ , "pro_significant"]] = "red"
    df$row_color_pro[is.na(df$name_pro)] = "grey"

    colnames(df) = gsub("rna_", "", gsub("pro_", "", colnames(df)))

    df = df[ , -grep("_vs_|significant", names(df))]
    dat = df[ , c(match("gene", colnames(df)), match("name_rna", colnames(df)), match("row_color_RNA", colnames(df)), 3:(match("name_rna", colnames(df)) -1), match("name_pro", colnames(df)), match("row_color_pro", colnames(df)), (match("name_rna", colnames(df)) +1):(match("name_pro", colnames(df)) -1))]
    names(dat) = gsub("\\.1","",colnames(dat))
  }

  if(if_chooseToshow) {
    if(length(chooseToshow) >= 1) {
      if(length(chooseToshow) == 1) {
        df1 = dat[match(chooseToshow, dat$gene), ]
        df1 = as.data.frame(df1)
        rownames(df1) = dat$gene[match(chooseToshow, dat$gene)]
        dat = df1
      } else {
        dat = dat[match(chooseToshow, dat$gene), ]
      }

    } else {
      return("at least one proteins should be selected")
    }
  }
  ######################################### next
  RNA = dat[ , match("row_color_RNA", colnames(dat)) : (match("name_pro", colnames(dat)) - 1)]
  PRO = dat[ , match("row_color_pro", colnames(dat)) : ncol(dat)]
  rownames(RNA) = dat$gene
  rownames(PRO) = dat$gene

  ################### PRO heatmap
  dat_pro_NA = PRO[apply(PRO[ , -1], 1, function(i){all(is.na(i))}), ]
  dat_pro = PRO[apply(PRO[ , -1], 1, function(i){!all(is.na(i))}), ]

  ht1=Heatmap(dat_pro[ , -1] - rowMeans(dat_pro[ , -1]), col = colorRamp2(seq(-Numcol_pro, Numcol_pro, (Numcol_pro/5)), rev(brewer.pal(11, color))),
              na_col = "grey96",row_names_gp = gpar(fontsize = row_font_size),#, col = row_name_color_pro),
              column_names_gp = gpar(fontsize = col_font_size),
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T,
              heatmap_legend_param = list(color_bar = "continuous",title_position = "topleft",legend_width = unit(5, "cm")),
              name = "Protein", cluster_columns = cluster_columns, row_km = 1, column_title = "Protein centered Enrichment (Log2)", column_title_gp = gpar(fontsize = column_title_size))

  rw = row_order(ht1)
  dat_pro = dat_pro[rw, ]
  dat_pro = rbind(dat_pro, dat_pro_NA)
  row_name_color_pro = dat_pro$row_color_pro
  dat_pro = dat_pro[ ,-1]
  dat_pro = dat_pro - rowMeans(dat_pro)

  ht1=Heatmap(dat_pro, col = colorRamp2(seq(-Numcol_pro, Numcol_pro, (Numcol_pro/5)), rev(brewer.pal(11, color))),
              na_col = "grey96",row_names_gp = gpar(fontsize = row_font_size, col = row_name_color_pro),
              column_names_gp = gpar(fontsize = col_font_size),
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T,
              heatmap_legend_param = list(color_bar = "continuous",title_position = "topleft",legend_width = unit(5, "cm")),
              name = "Protein", cluster_columns = cluster_columns, row_km = 1, column_title = "Protein centered Enrichment (Log2)", column_title_gp = gpar(fontsize = column_title_size), cluster_rows = FALSE)

  ht_pro <- draw(ht1, heatmap_legend_side = c("right"))
  ################### RNA heatmap
  dat_rna = RNA[match(rownames(dat_pro), rownames(RNA)), ]
  row_name_color_rna = dat_rna$row_color_RNA
  dat_rna = dat_rna[ ,-1]
  dat_rna = dat_rna - rowMeans(dat_rna)


  ht2=Heatmap(dat_rna, col = colorRamp2(seq(-Numcol_rna, Numcol_rna, (Numcol_rna/5)), rev(brewer.pal(11, color))),
              na_col = "grey96",row_names_gp = gpar(fontsize = row_font_size, col = row_name_color_rna),
              column_names_gp = gpar(fontsize = col_font_size),
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T,
              heatmap_legend_param = list(color_bar = "continuous",title_position = "topleft"), name = "RNA", cluster_columns = cluster_columns, row_km = 1,cluster_rows = F, column_title = "RNA centered normalized count", column_title_gp = gpar(fontsize = column_title_size))#width = unit(6, "cm")

  ht_rna <- draw(ht2, heatmap_legend_side = c("right"))

  ht_lis = ht1 + ht2
  hhh = draw(ht_lis, ht_gap = unit(0.4, "cm"), main_heatmap = "Protein", auto_adjust = FALSE,legend_title_position = "top",heatmap_legend_side = "right" )

    return(lis = list(ht_pro = ht_pro, ht_rna = ht_rna, ht_merge = hhh, dat = dat))
  }


# ## speed up the speed of GO analysis
# kegg_list <- function (db) {
#   url <- paste0("http://rest.kegg.jp/list/", db, collapse = "")
#   clusterProfiler:::kegg_rest(url)
# }
#
# ko2name <- function (ko) {
#   p <- clusterProfiler:::kegg_list("pathway")
#   ko2 <- gsub("^ko", "path:map", ko)
#   ko.df <- data.frame(ko = ko, from = ko2)
#   res <- merge(ko.df, p, by = "from", all.x = TRUE)
#   res <- res[, c("ko", "to")]
#   colnames(res) <- c("ko", "name")
#   return(res)
# }
#
#
# get_GOdata_Env <- function () {
#   if (!exists(".GO_data_Env", envir = .GlobalEnv)) {
#     pos <- 1
#     envir <- as.environment(pos)
#     assign(".GO_data_Env", new.env(), envir=envir)
#   }
#   get(".GO_data_Env", envir = .GlobalEnv)
# }
#
# enrichGO <- function(gene,
#                      OrgDb,
#                      keyType = "ENTREZID",
#                      ont="MF",
#                      pvalueCutoff=0.05,
#                      pAdjustMethod="BH",
#                      universe,
#                      qvalueCutoff = 0.2,
#                      minGSSize = 10,
#                      maxGSSize = 500,
#                      readable=FALSE, pool=FALSE) {
#   ont %<>% toupper
#   ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))
#
#   OrgDb_name = OrgDb$packageName
#   OrgDb_version = package.version(OrgDb$packageName)
#   GO_DATAfolder = paste(system.file(package = "DEP2"),"/GOdata",sep="")
#   if (!file.exists(GO_DATAfolder)){
#     dir.create(file.path(GO_DATAfolder))
#   }
#   exist_GOdata =  list.files(file.path(GO_DATAfolder))
#   GO_DATAfile = paste(OrgDb_name,"_",OrgDb_version,"_",keyType,".RDS",sep = "")
#
#   GO_DATA_environment = get_GOdata_Env()
#
#   if(GO_DATAfile %in% exist_GOdata){
#     use_cashed = F  ## load a exiting RDS
#     if(exists("GO_DATAfile",envir = GO_DATA_environment) & exists("GO_DATA",envir = GO_DATA_environment)){
#       GO_DATAfile2 = get("GO_DATAfile",envir = GO_DATA_environment)
#       if(GO_DATAfile2 == GO_DATAfile)
#         use_cashed = T  ## a correct RDS was loaded
#       # GO_DATA = get("GO_DATA",envir = GO_DATA_environment)
#     }
#
#     if(use_cashed){
#       GO_DATA = get("GO_DATA",envir = GO_DATA_environment)
#     }else{
#       GO_DATA = readRDS(paste(GO_DATAfolder,"/",GO_DATAfile,sep = ""))
#       assign("GO_DATA",GO_DATA,envir = GO_DATA_environment)
#       assign("GO_DATAfile",GO_DATAfile,envir = GO_DATA_environment)
#
#       cat("loaded a existing GOdata\n")
#     }
#
#   }else{
#     ## creat a GO_DATA for particular GO options: OrgDb,ont,keyType
#     GO_DATA <- list(ALL = as.list(get_GO_data(OrgDb, "ALL", keyType)),
#                     BP = as.list(get_GO_data(OrgDb, "BP", keyType)),
#                     MF = as.list(get_GO_data(OrgDb, "MF", keyType)),
#                     CC = as.list(get_GO_data(OrgDb, "CC", keyType)))
#     ## save the GO_DATA as an rds file
#     saveRDS(GO_DATA,file = paste(GO_DATAfolder,"/",GO_DATAfile,sep = ""))
#
#     assign("GO_DATA",GO_DATA,envir = GO_DATA_environment)
#     assign("GO_DATAfile",GO_DATAfile,envir = GO_DATA_environment)
#
#     cat(paste("saved a new GOdata"))
#   }
#
#   GO_DATA = as.environment(GO_DATA[[ont]])
#   # cat(ont)
#   # cat("GO_DATA finished\n\n")
#   # GO_DATA <- get_GO_data(OrgDb, ont, keyType)
#
#   if (missing(universe))
#     universe <- NULL
#
#   if (ont == "ALL" && !pool) {
#     lres <- lapply(c("BP", "CC", "MF"), function(ont)
#       suppressMessages(enrichGO(gene, OrgDb, keyType, ont,
#                                 pvalueCutoff, pAdjustMethod, universe,
#                                 qvalueCutoff, minGSSize, maxGSSize
#       ))
#     )
#
#     lres <- lres[!vapply(lres, is.null, logical(1))]
#     if (length(lres) == 0)
#       return(NULL)
#
#     df <- do.call('rbind', lapply(lres, as.data.frame))
#     geneSets <- lres[[1]]@geneSets
#     if (length(lres) > 1) {
#       for (i in 2:length(lres)) {
#         geneSets <- append(geneSets, lres[[i]]@geneSets)
#       }
#     }
#     res <- lres[[1]]
#     res@result <- df
#     res@geneSets <- geneSets
#   } else {
#     res <- enricher_internal(gene,
#                              pvalueCutoff=pvalueCutoff,
#                              pAdjustMethod=pAdjustMethod,
#                              universe = universe,
#                              qvalueCutoff = qvalueCutoff,
#                              minGSSize = minGSSize,
#                              maxGSSize = maxGSSize,
#                              USER_DATA = GO_DATA
#     )
#
#     if (is.null(res))
#       return(res)
#   }
#
#   res@keytype <- keyType
#   res@organism <- get_organism(OrgDb)
#   if(readable) {
#     res <- setReadable(res, OrgDb)
#   }
#   res@ontology <- ont
#
#   if (ont == "ALL") {
#     res <- add_GO_Ontology(res, GO_DATA)
#   }
#   return(res)
# }
#
# environment(enrichGO) = asNamespace("clusterProfiler")
# assignInNamespace("enrichGO", value = enrichGO, ns = "clusterProfiler")



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

## load PPIdata of spcies from DEP/PPIdata dir to strEnv=.string_Env
load_PPIdata = function(speciesname){
  strEnv <- get_string_Env()
  if(exists(c("load succed"),envir = strEnv) && get("speciesname",envir = strEnv) == speciesname && get("load succed",envir = strEnv)){
    return("use_cache")
  }

  scriptdir <- system.file("PPIdata",package = "DEP2")
  version <- "11.0"

  speciesID = switch (speciesname,
                      "Anopheles" = "7165",
                      "Arabidopsis" = "3702",
                      "Bovine" = "9913",
                      "Worm" = "6239",
                      "Canine" = "9612",
                      "Fly" = "7227",
                      "Zebrafish" = "7955",
                      "E coli strain K12" = "511145",
                      "Chicken" = "9031",
                      "Human" = "9606",
                      "Mouse" = "10090",
                      "Rhesus" = "9544",
                      "Malaria" = "5833",
                      "Chimp" = "9598",
                      "Rat" = "10116",
                      "Yeast" = "4932",
                      "Streptomyces coelicolor" = "100226",
                      "Pig" = "9823",
                      "Toxoplasma gondii" = "5811",
                      "Xenopus" = "8364"
  )

  cat(paste("PPIdata directory is",scriptdir,"\n"))
  cat(paste("Your choice is ",speciesname ,"\t", "The species ID is ",speciesID,"\n",sep=""))

  ###加载所需文件到string_Env中
  proteinsDffile <- paste(scriptdir,"/",speciesname, "/", speciesID, ".protein.info.v", version, ".txt.gz", sep="")
  aliasesfile = paste(scriptdir,"/",speciesname, "/", speciesID, ".protein.aliases.v", version, ".txt.gz", sep="")
  links_detailfile = paste(scriptdir,"/",speciesname, "/", speciesID, ".protein.links.detailed.v", version, ".txt.gz", sep="")

  # existfiles = list.files(file.path(scriptdir,speciesname))

  if(!all( file.exists(c(proteinsDffile,aliasesfile,links_detailfile)) )){
    stop(paste("Lack the necessary String datasets in ",file.path(scriptdir,speciesname),
               "\n Please download the info,aliases,detailed txt.gz files of this species from String, and paste to ",
               file.path(scriptdir,speciesname),sep = ""))
  }


  assign("proteinsDf",fread(proteinsDffile, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE, quote=""),envir = strEnv)
  assign("aliasDf",suppressWarnings(fread(aliasesfile)),envir = strEnv)
  assign("protein_links_detail",fread(links_detailfile),envir = strEnv)
  assign("speciesID",speciesID,envir= strEnv)
  assign("speciesname",speciesname,envir= strEnv)
  assign("load succed",T,envir=strEnv)
  cat("load succeed")

}

myget_interactions = function(protein_detail,mapped_data,choose_scores,score_cutoff = 400){
  protein_detail2 <- as.data.frame(protein_detail)[,c("protein1","protein2",choose_scores)]
  protein_detail2 <- dplyr::filter(protein_detail2,rowMax(as.matrix(protein_detail2[,choose_scores]))>=score_cutoff)
  if(!"combined_score" %in% choose_scores){
    mat <- as.matrix(protein_detail2[,choose_scores])
    mat2 <- (1000-mat)/1000
    # library(parallel)
    # cl <- makeCluster(3)
    # protein_detail2$custom_combine_score <- 1000 - 1000*parApply(cl=cl,mat2,1,prod)
    # stopCluster(cl)
    # custom_combine_score <- 1-s
    protein_detail2$combined_score <- 1000 - 1000*apply(mat2,1,prod)
  }

  string_ids <- mapped_data$STRING_id
  link_table <- protein_detail2 %>% dplyr::filter((protein1 %in% string_ids) & (protein2 %in% string_ids))
  colnames(link_table)[1:2] <- c("from","to")
  return(link_table)
}


get_proteins = function(proteinsDf){
  '
Description:
  Returns the STRING proteins data frame.
  (it downloads and caches the information the first time that is called).

Author(s):
   Andrea Franceschini
'
  # temp = downloadAbsentFile(paste("https://stringdb-static.org/download/protein.info.v", version, "/", species, ".protein.info.v", version, ".txt.gz", sep=""), oD=input_directory)
  # temp <- proteinsDf
  # proteinsDf <- fread(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE, quote="")
  proteinsDf2 = subset(proteinsDf, select=c("protein_external_id",  "preferred_name", "protein_size", "annotation"))
  proteins <- proteinsDf2
  # cat(paste("load proteinsDFfile: ",temp ,"\n", sep=""))

  return(proteins)
}

get_aliases = function(aliasDf,proteinsDf){
  '
Description:
  Loads and returns STRING alias table.

Author(s):
   Andrea Franceschini
'
  # temp = aliasesfile

  proteins <- get_proteins(proteinsDf)
  colnames(aliasDf) <- c("STRING_id", "alias", "sources")
  aliasDf = subset(aliasDf, select=c("STRING_id", "alias"))
  pr1=data.frame(STRING_id=proteins$protein_external_id, alias=proteins$preferred_name, stringsAsFactors=FALSE)
  pr2=data.frame(STRING_id=proteins$protein_external_id, alias=proteins$protein_external_id, stringsAsFactors=FALSE)
  pr3=data.frame(STRING_id=proteins$protein_external_id, alias=unlist(strsplit(proteins$protein_external_id, "\\."))[seq(from=2, to=2*nrow(proteins), by=2)], stringsAsFactors=FALSE)
  #if(takeFirst){aliasDf = subset(aliasDf, !(alias %in% proteins$preferred_name) & !(alias %in% proteins$protein_external_id) )  }
  aliasDf = subset(aliasDf, !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(proteins$preferred_name)) &
                     !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(proteins$protein_external_id))  &
                     !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(unlist(strsplit(proteins$protein_external_id, "\\."))[seq(from=2, to=2*nrow(proteins), by=2)])) )

  aliasDf2=rbind(pr1,pr2,pr3, aliasDf)
  aliases_tf <<- aliasDf2

  return(aliasDf2)
}

merge.with.order <- function(x,y, ..., sort = T)
{
  # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
  order.by.id...and.remove.it <- function(DATA)
  {
    # gets in a data.frame with the "id..." column.  Orders by it and returns it
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")

    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }

  if(sort==F){ return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
  } else {return(merge(x=x,y=y,..., sort = sort))}

}

multi_map_df <- function(dfToMap, dfMap, strColsFrom, strColFromDfMap, strColToDfMap, caseSensitive=FALSE){

  tempMatr = matrix(NA, length(strColsFrom), nrow(dfToMap))
  for(i in 1:length(strColsFrom)){
    if(!caseSensitive){
      tempMatr[i,] = as.vector(dfToMap[,strColsFrom[i]])
      dfToMap[,strColsFrom[i]] = toupper(iconv(dfToMap[,strColsFrom[i]], "WINDOWS-1252", "UTF-8"))
    }
  }
  #error
  dfMap = as.data.frame(dfMap)
  if(!caseSensitive){ dfMap[,strColFromDfMap] = toupper(iconv(dfMap[,strColFromDfMap], "WINDOWS-1252", "UTF-8")) }

  dfMap2 = unique(subset(dfMap, select=c(strColFromDfMap, strColToDfMap)))
  df2 = merge.with.order(dfToMap, dfMap2, by.x=strColsFrom[1], by.y=strColFromDfMap, all.x=TRUE, sort=FALSE)
  if(length(strColsFrom) > 1){
    for(i in 2:length(strColsFrom)){
      dfna = delColDf(subset(df2, is.na(as.vector(df2[, strColToDfMap]))), strColToDfMap)
      dfgood = subset(df2, !is.na(as.vector(df2[, strColToDfMap])))
      df3 = merge.with.order(dfna, dfMap2, by.x=strColsFrom[i], by.y=strColFromDfMap, all.x=TRUE, sort=FALSE)
      df2 = rbind(dfgood, df3)
    }
  }

  for(i in 1:length(strColsFrom)){
    if(!caseSensitive && length(tempMatr[i,])==length(df2[,strColsFrom[i]])) df2[,strColsFrom[i]] = tempMatr[i,]
  }

  return(df2)
}

mymap = function(my_data_frame,
                 my_data_frame_id_col_names,proteinsDf,aliasDf,
                 takeFirst=TRUE, removeUnmappedRows=T, quiet=FALSE,reload=F
){
  '
Description:
  Maps the gene identifiers of the input dataframe to STRING identifiers.
  It returns the input dataframe with the "STRING_id" additional column.

Input parameters:
  "my_data_frame"                 data frame provided as input.
  "my_data_frame_id_col_names"    vector contatining the names of the columns of "my_data_frame" that have to be used for the mapping.
  "takeFirst"                     boolean indicating what to do in case of multiple STRING proteins that map to the same name.
                                      If TRUE, only the first of those is taken. Otherwise all of them are used. (default TRUE)
  "removeUnmappedRows"            remove the rows that cannot be mapped to STRING
                                      (by default those lines are left and their STRING_id is set to NA)
  "quiet"                         Setting this variable to TRUE we can avoid printing the warning relative to the unmapped values.

Author(s):
   Andrea Franceschini
'
  aliasDf2=get_aliases(aliasDf=aliasDf,proteinsDf=proteinsDf)

  tempDf = multi_map_df(my_data_frame, aliasDf2, my_data_frame_id_col_names, "alias", "STRING_id")
  naDf = subset(tempDf, is.na(STRING_id))
  if(nrow(naDf) > 0 & !quiet) cat(paste("We couldn't map to STRING ", as.integer((nrow(naDf)/nrow(tempDf))*100), "% of your identifiers" , sep=""))
  if(removeUnmappedRows) tempDf = subset(tempDf, !is.na(STRING_id))

  return(tempDf)
}




filterlink <- function(linkstable,scorecutoff){
  linkstable2 <- linkstable %>% dplyr::filter(combined_score>=scorecutoff)
  for (i in 1:(nrow(linkstable2)-1)) {
    if(i >= nrow(linkstable2) ) break
    for (j in (i+1):nrow(linkstable2)) {
      if(linkstable2[i,1]==linkstable2[j,2] &&linkstable2[i,2]==linkstable2[j,1]){
        linkstable2 = linkstable2[-j,]
        break
      }
    }
  }
  return(linkstable2)
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
  links4 <<- links3;nodes4 <<-nodes
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


get_dep_siglist = function (dep, alpha = 0.05, lfc = 1, threshold_method = c("intersect","curve"), curvature=0.6,x0_fold = 2)
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

  genelist_env = get_dep_genelist_Env()
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
        assign(value = temp, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Sig",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)
      } else if(threshold_method == "curve") {
        assign(value = temp, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Sig",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Sig",sep = "_"),.), sep = "") ,
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
        assign(value = temp_up, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Up",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_up, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Up",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)

      } else if(threshold_method == "curve") {
        assign(value = temp_up, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Up",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_up, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Up",sep = "_"),.), sep = "") ,
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
        assign(value = temp_dn, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_dn, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,alpha,lfc,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)
      } else if(threshold_method == "curve") {
        assign(value = temp_dn, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_env)
        assign(value = temp_dn, x = paste("DEP_",colnames(res)[cols_p[i]] %>% gsub("p.adj",paste(threshold_method,curvature,x0_fold,"Dn",sep = "_"),.), sep = "") ,
               envir = genelist_all_env)
      }

    }
  }
  return(ls(envir = get_dep_genelist_Env()))

}


generate_dep_siglist = function(dep, alpha = 0.05, lfc = 1,compare,type=c("Up","Dn","Sig"), threshold_method = c("intersect","curve"), curvature=0.6,x0_fold = 2){
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

  if(length(grep(paste(compare,"_p.adj",sep = ""),colnames(res))) ==0){
    cat("compare do not exit")
  }

  if( length(threshold_method)!=1 | (!threshold_method %in% c("intersect","curve")) ){
    stop("threshold_method should be one of 'intersect'/'curve'",call. = F)
  }

  if(threshold_method == "intersect") {
    if(type == "Up"){
      tem = (res[,paste(compare,"_p.adj",sep = "")]<alpha) & (res[,paste(compare,"_ratio",sep="")] > lfc)
    }else if(type == "Dn"){
      tem = (res[,paste(compare,"_p.adj",sep = "")]<alpha) & (res[,paste(compare,"_ratio",sep="")] < -lfc)
    }else if(type == "Sig"){
      tem = (res[,paste(compare,"_p.adj",sep = "")]<alpha) & (abs(res[,paste(compare,"_ratio",sep="")]) > lfc)
    }
  }

  if(threshold_method == "curve") {
    x = res[, paste(compare,"_ratio",sep="")]
    x = fun.outlier(x) %>% na.omit()
    fit <- fitnormal(x)
    σ = sqrt(fit$theta[2])
    x0 = x0_fold*σ
    polar = ifelse((res[,paste(compare,"_ratio",sep="")] > 0),1,-1)
    if(type == "Up") {
      tem <- ( -log10(res[,paste(compare,"_p.adj",sep = "")]) > curvature/abs(res[,paste(compare,"_ratio",sep="")] - polar*x0) &
          (res[,paste(compare,"_ratio",sep="")] >= 0) & res[,paste(compare,"_ratio",sep="")]>x0 )

    } else if(type == "Dn") {
      tem <- ( -log10(res[,paste(compare,"_p.adj",sep = "")]) > curvature/abs(res[,paste(compare,"_ratio",sep="")] - polar*x0) &
                 (res[,paste(compare,"_ratio",sep="")] < 0) & res[,paste(compare,"_ratio",sep="")] < -x0 )

    } else if(type == "Sig") {
      tem <- ( -log10(res[,paste(compare,"_p.adj",sep = "")]) > curvature/abs(res[,paste(compare,"_ratio",sep="")] - polar*x0) &
                 ifelse((res[,paste(compare,"_ratio",sep="")] >= 0),res[,paste(compare,"_ratio",sep="")]>x0,res[,paste(compare,"_ratio",sep="")] < -x0) )

    }
  }

  # the significant column is according to shiny global variables p.adj and lfc cutoff, not the p.adj and lfc cutoff of generate genelist
  temp = res[which(tem),c("name","ID",
                          paste(compare,"_ratio",sep=""),paste(compare,"_p.adj",sep = ""),
                          paste(compare,"_significant",sep = ""))]


  if(nrow(temp)>0){
    colnames(temp) = c("symbol","ID", "L2FC", "p.adj","significant")
    temp$significant = TRUE
    temp = temp[ , c(1, 3, 2, 4, 5)]
    # assign(value = temp, x = paste("DEP",compare,alpha,lfc,type, sep = "_") ,
    #        envir = genelist_env)
    return(temp)
  }else{return(NULL)}
}

MyVenn <- function(x){
  venn <- Venn(x)
  data <- process_data(venn)
  items <- data@region %>% dplyr::rowwise() %>%
    dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, collapse = " "), width = 40)) %>%
    sf::st_as_sf() %>% dplyr::mutate(ratio = round(count/sum(data@region$count),3)) %>% dplyr::mutate(count2 = paste(count,"(",ratio*100,"%)",sep = ""))
  label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
  p <- ggplot(items) + geom_sf(aes_string(fill = "count"),lwd=0.5,color ="grey80") +
    geom_sf_text(aes_string(label = "name"), data = data@setLabel,inherit.aes = F) +
    geom_text(aes_string(label = "count2", text = "text"), x = label_coord[, 1], y = label_coord[,2], show.legend = FALSE,size=3) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme_void()
  p
  ax <- list(showline = FALSE)
  p <- plotly::ggplotly(p, tooltip = c("text")) %>% plotly::layout(xaxis = ax,yaxis = ax)
  return(list(p = p, df = items %>% as.data.frame() %>% dplyr::select(., name, count, ratio, count2, text) %>% dplyr::rename(., count_ratio = count2) ))
}

# get dep all list of all contrasts of dep,1:2columns gene and log2 fc
get_dep_listforgsea = function (dep)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
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
  cols_p <- grep("_p.adj", colnames(res))
  cols_diff <- grep("_ratio", colnames(res))
  cols_sig <- grep("_significant", colnames(res))


  genelist_gsea_env = get_gsea_genelist_Env()


  for(i in 1:length(cols_p)){
    temp = res[,c(1:2,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp)>0){
      colnames(temp) = c("symbol","ID","L2FC","p.adj","significant")
      temp = temp[ , c(1, 3, 2, 4, 5)]
      assign(value = temp, x = paste("DEP_",colnames(res)[cols_p[i]], sep = "") %>% gsub("_p.adj", "", .) ,
             envir = genelist_gsea_env)
    }

  }
  return(ls(envir = get_gsea_genelist_Env()))

}

### get all significant and up, down result of contrasts for DEG-RNAseq, and assign to deg evir, used for genelist tools to refresh the exited gene list and the source drag box of compare
# deg:  data frame, the result of get_res and add_rejection for DEG-RNAseq options. set test type to all is the best, in order to provide more choices for user
# alpha: the p adj cutoff of DEG-RNAseq options
# lfc: the log2 fold change cutoff of DEG-RNAseq options
get_deg_siglist = function (deg, alpha = 0.05, lfc = 1)
{
  if (is.integer(alpha))
    alpha <- as.numeric(alpha)
  if (is.integer(lfc))
    lfc <- as.numeric(lfc)
  assertthat::assert_that(is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                          length(lfc) == 1)
  # row_data <- rowData(dep, use.names = FALSE) %>% as.data.frame()
  res = deg %>% as.data.frame(.) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname)
  if (any(!c("Gene") %in% colnames(res))) {
    stop("'Gene' column are not present in '",
         deparse(substitute(deg)), "'\nRun get_res() to obtain the required column",
         call. = FALSE)
  }
  if (length(grep("_padj|_log2FoldChange", colnames(res))) < 1) {
    stop("'[contrast]_log2FoldChange' and/or '[contrast]_padj' columns are not present in '",
         deparse(substitute(deg)), "'\nRun get_res() to obtain the required columns",
         call. = FALSE)
  }
  cols_p <- grep("_padj", colnames(res))
  cols_diff <- grep("_log2FoldChange", colnames(res))
  cols_sig <- grep("_significant", colnames(res))

  genelist_env = get_deg_genelist_Env()
  genelist_all_env = get_all_genelist_Env()
  rm(list = ls(envir = genelist_env),envir = genelist_env)
  # rm(list = ls(envir = genelist_all_env),envir = genelist_all_env)

  for(i in 1:length(cols_p)){
    tem = res[,cols_p[i]]<alpha & abs(res[,cols_diff[i]] )> lfc
    ind = which(colnames(res) %in% c("Gene", "symbol"))
    temp = res[which(tem),c(ind,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp)>0){
      if(length(ind) == 1) {
        colnames(temp) = c("Gene","L2FC","p.adj","significant")
        temp$significant = TRUE
        assign(value = temp, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Sig",sep = "_"),.), sep = "") , envir = genelist_env)
        assign(value = temp, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Sig",sep = "_"),.), sep = "") , envir = genelist_all_env)
      }
      if(length(ind) == 2) {
        colnames(temp) = c("Gene","symbol","L2FC","p.adj","significant")
        temp$significant = TRUE
        temp = temp[ , c("symbol", "L2FC", "Gene","p.adj","significant")]
        assign(value = temp, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Sig",sep = "_"),.), sep = "") , envir = genelist_env)
        assign(value = temp, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Sig",sep = "_"),.), sep = "") , envir = genelist_all_env)

      }
    }

    tem_up = (res[,cols_p[i]]<alpha) & (res[,cols_diff[i]] > lfc)
    temp_up = res[which(tem_up),c(ind,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp_up)>0){
      if(length(ind) == 1) {
        colnames(temp_up) = c("Gene","L2FC","p.adj","significant")
        temp_up$significant = TRUE
        assign(value = temp_up, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Up",sep = "_"),.), sep = "") , envir = genelist_env)
        assign(value = temp_up, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Up",sep = "_"),.), sep = "") , envir = genelist_all_env)
      }
      if(length(ind) == 2) {
        colnames(temp_up) = c("Gene","symbol","L2FC","p.adj","significant")
        temp_up$significant = TRUE
        temp_up = temp_up[ , c("symbol", "L2FC", "Gene", "p.adj","significant")]
        assign(value = temp_up, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Up",sep = "_"),.), sep = "") , envir = genelist_env)
        assign(value = temp_up, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Up",sep = "_"),.), sep = "") , envir = genelist_all_env)

      }
    }

    tem_dn = (res[,cols_p[i]]<alpha) & (res[,cols_diff[i]] < -lfc)
    temp_dn = res[which(tem_dn),c(ind,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp_dn)>0){
      if(length(ind) == 1) {
        colnames(temp_dn) = c("Gene","L2FC","p.adj","significant")
        temp_dn$significant = TRUE
        assign(value = temp_dn, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Dn",sep = "_"),.), sep = "") , envir = genelist_env)
        assign(value = temp_dn, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Dn",sep = "_"),.), sep = "") , envir = genelist_all_env)
      }
      if(length(ind) == 2) {
        colnames(temp_dn) = c("Gene","symbol","L2FC","p.adj","significant")
        temp_dn$significant = TRUE
        temp_dn = temp_dn[ , c("symbol", "L2FC", "Gene", "p.adj","significant")]
        assign(value = temp_dn, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Dn",sep = "_"),.), sep = "") , envir = genelist_env)
        assign(value = temp_dn, x = paste("DEG_",colnames(res)[cols_p[i]] %>% gsub("padj",paste(alpha,lfc,"Dn",sep = "_"),.), sep = "") , envir = genelist_all_env)

      }
    }
  }
  return(ls(envir = get_deg_genelist_Env()))

}



##############################
### get significant, up or down by custom p.adj and lfc cutoff result of a contrast for DEG-RNAseq, it is used for genelist tools to generate a gene list , and then are assigned to deg evir
# deg:  data frame, the result of get_res and add_rejection for DEG-RNAseq options. set test type to all is the best, in order to provide more choices for user
# alpha: the p adj cutoff of DEG-RNAseq options in the generate genelist by cutoff panel of genelist tool
# lfc: the log2 fold change cutoff of DEG-RNAseq options in the generate genelist by cutoff panel of genelist tool
# compare: the contrast to generate gene list, the choices are all options in deg
# type: the type to generate gene list. one of up, Fn, and sig. Up: up regulated, Dn: down regulated, sig: up and down regulated
generate_deg_siglist = function(deg, alpha = 0.05, lfc = 1,compare,type=c("Up","Dn","Sig")){
  if (is.integer(alpha))
    alpha <- as.numeric(alpha)
  if (is.integer(lfc))
    lfc <- as.numeric(lfc)
  assertthat::assert_that(is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                          length(lfc) == 1)
  # row_data <- rowData(dep, use.names = FALSE) %>% as.data.frame()
  res = deg %>% as.data.frame(.) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname)
  if (any(!c("Gene") %in% colnames(res))) {
    stop("'Gene' column are not present in '",
         deparse(substitute(deg)), "'\nRun get_res() to obtain the required column",
         call. = FALSE)
  }
  if (length(grep("_padj|_log2FoldChange", colnames(res))) < 1) {
    stop("'[contrast]_log2FoldChange' and/or '[contrast]_padj' columns are not present in '",
         deparse(substitute(deg)), "'\nRun get_res() to obtain the required columns",
         call. = FALSE)
  }

  if(length(grep(paste(compare,"_padj",sep = ""),colnames(res))) ==0){
    cat("compare do not exit")
  }
  if(type == "Up"){
    tem = (res[,paste(compare,"_padj",sep = "")]<alpha) & (res[,paste(compare,"_log2FoldChange",sep="")] > lfc)
  }else if(type == "Dn"){
    tem = (res[,paste(compare,"_padj",sep = "")]<alpha) & (res[,paste(compare,"_log2FoldChange",sep="")] < -lfc)
  }else if(type == "Sig"){
    tem = (res[,paste(compare,"_padj",sep = "")]<alpha) & (abs(res[,paste(compare,"_log2FoldChange",sep="")]) > lfc)
  }

  tem[which(is.na(tem))] = FALSE

  ind =  colnames(res)[which(colnames(res) %in% c("Gene", "symbol"))]

  temp = res[which(tem),c(ind,
                   paste(compare,"_log2FoldChange",sep=""),paste(compare,"_padj",sep = ""),
                   paste(compare,"_significant",sep = ""))]

  if(nrow(temp)>0){
    if(length(ind) == 1) {
      colnames(temp) = c("Gene", "L2FC", "p.adj","significant")
    }

    if(length(ind) == 2) {
      colnames(temp) = c("Gene","symbol", "L2FC", "p.adj","significant")
      temp = temp[ , c("symbol", "L2FC", "Gene", "p.adj","significant")]
    }

    temp$significant = TRUE
    # assign(value = temp, x = paste("DEG",compare,alpha,lfc,type, sep = "_") ,
    #        envir = genelist_env)
    return(temp)
  }else{return(NULL)}
}

# get deg(DEG-RNAseq) all list of all contrasts,1:2columns gene and log2 fc
get_deg_listforgsea = function (deg)
{
  res = deg %>% as.data.frame(.) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname)
  if (any(!c("Gene") %in% colnames(res))) {
    stop("'Gene' column are not present in '",
         deparse(substitute(deg)), "'\nRun get_res() to obtain the required column",
         call. = FALSE)
  }
  if (length(grep("_padj|_log2FoldChange", colnames(res))) < 1) {
    stop("'[contrast]_log2FoldChange' and/or '[contrast]_padj' columns are not present in '",
         deparse(substitute(deg)), "'\nRun get_res() to obtain the required columns",
         call. = FALSE)
  }
  cols_p <- grep("_padj", colnames(res))
  cols_diff <- grep("_log2FoldChange", colnames(res))
  cols_sig <- grep("_significant", colnames(res))

  genelist_gsea_env = get_gsea_genelist_Env()


  for(i in 1:length(cols_p)){
    ind = which(colnames(res) %in% c("Gene", "symbol"))
    temp = res[,c(ind,cols_diff[i],cols_p[i],cols_sig[i])]
    if(nrow(temp)>0){
      if(length(ind) == 1) {
        colnames(temp) = c("Gene","L2FC","p.adj","significant")
        temp = temp[which(!is.na(temp[,1])), ]
        assign(value = temp, x = paste("DEG_",colnames(res)[cols_p[i]], sep = "") %>% gsub("_padj", "", .) , envir = genelist_gsea_env)
      }
      if(length(ind) == 2) {
        colnames(temp) = c("Gene","symbol","L2FC","p.adj","significant")
        temp = temp[ , c("symbol", "L2FC", "Gene","p.adj","significant")]
        temp = temp[which(!is.na(temp[,1])), ]
        assign(value = temp, x = paste("DEG_",colnames(res)[cols_p[i]], sep = "") %>% gsub("_padj", "", .) , envir = genelist_gsea_env)
      }
    }
  }
  return(ls(envir = get_gsea_genelist_Env()))

}

get_all_lists <- function(){
  c(ls(envir = get_dep_genelist_Env()),
    ls(envir = get_deg_genelist_Env()),
    ls(envir = get_imported_genelist_Env()))
}


# get_all_lists <- function(){
#   ls(envir = get_all_genelist_Env())
# }

get_all_lists_for_gsea <- function(){
  ls(envir = get_gsea_genelist_Env())
}

# get_mergelists = function(listnames,env_list){
#   lists = env_list %>% sapply(.,function(x){ls(envir = x)})
#   names(lists) = paste(1: length(lists),"_",sep = "")
#   lists = unlist(lists)
#   names(lists) = names(lists) %>% sapply(., function(x){strsplit(x,"_")[[1]][1]})

#   mergelist = vector()
#   for(i in listnames){
#     theenv =  names(lists)[match(i,lists)][1] %>% as.numeric()
#     theenv2 = env_list[[ names(lists)[match(i,lists)][1] %>% as.numeric() ]]
#     thelist = get(i,envir = theenv2)[,1]
#     if(names(lists)[match(i,lists)][1] %>% as.numeric() == 1){
#       thelist = thelist %>% sapply(., function(x){strsplit(x,"\\.")[[1]][1]}) %>% unique
#     }
#     mergelist = c(mergelist,thelist) %>% unique()
#   }
#   return(mergelist)
# }


# for annotation, ORA, GSEA, PPI, DEG-RNAseq select organism
# annoSpecies_df <<-
#     data.frame(
#       species = c(
#         "", "Anopheles", "Arabidopsis", "Bovine", "Worm",
#         "Canine", "Fly", "Zebrafish", "E coli strain K12",
#         "E coli strain Sakai", "Chicken", "Human", "Mouse",
#         "Rhesus", "Malaria", "Chimp", "Rat",
#         "Yeast", "Streptomyces coelicolor", "Pig", "Toxoplasma gondii",
#         "Xenopus"
#       ),
#       pkg = c(
#         "", "org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
#         "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
#         "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
#         "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
#         "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db",
#         "org.Xl.eg.db"
#       ),
#       stringsAsFactors = FALSE
#     )
# annoSpecies_df$organism = c("", "aga", "ath", "bta", "cel", "cfa", "dme", "dre", "eco", "ecs", "gga", "hsa", "mmu", "mcc", "pfa", "ptr", "rno", "sce", "sco", "ssc", "tgo", "xla")
# annoSpecies_df_for_reactome <- annoSpecies_df[c(12,17,13,5,18,8,7), ]
# annoSpecies_df_for_reactome$organism = c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly")
# annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species), ]
# rownames(annoSpecies_df) <- annoSpecies_df$species # easier to access afterwards




################### for DEP-LFQ significant cutoff method curve and intersect
## rm outlier
fun.outlier <- function(x,time.iqr=1.5) {
  outlier.low <- quantile(x,probs=c(0.25))-IQR(x)*time.iqr
  outlier.high <- quantile(x,probs=c(0.75))+IQR(x)*time.iqr
  x[which(x>outlier.high | x<outlier.low)]<-NA
  x
}

## fitnormal function
fitnormal = function (x, exact = TRUE) {
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




#' Mark significant proteins
#'
#' \code{add_rejections} marks significant proteins based on defined cutoffs.
#'
#' @param diff SummarizedExperiment,
#' Proteomics dataset on which differential enrichment analysis
#' has been performed (output from \code{\link{test_diff}()}).
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
#' @export
add_rejections <- function (diff, alpha = 0.05, lfc = 1,thresholdmethod="intersect",curvature=1,x0_fold = 2)
{

  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(diff)), "'\nRun make_unique() and make_se() to obtain the required columns",
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
      x = fun.outlier(x) %>% na.omit()
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
#################

########################################################################### For DEP-LFQ part, add two plots after pca plot
library(umap)
library(Rtsne)
#indicate can be one , two, three element, the third means facet_wrap
# if_square: if xlim == ylim, default FALSE
plot_umap <- function (dep, indicate = c("condition",
                                         "replicate"), label = FALSE, n = nrow(dep), point_size = 4,
                       label_size = 3, plot = TRUE, if_square = FALSE,
                       n_neighbors=4,random_state=10086)
{
  # if (is.integer(x))
  #   x <- as.numeric(x)
  # if (is.integer(y))
  #   y <- as.numeric(y)
  if (is.integer(n))
    n <- as.numeric(n)
  if (is.integer(point_size))
    point_size <- as.numeric(point_size)
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.numeric(n), length(n) == 1, is.character(indicate),
                          is.logical(label), is.numeric(point_size), length(point_size) ==1,
                          is.numeric(label_size), length(label_size) ==1,
                          is.logical(plot), length(plot) == 1,
                          is.numeric(n_neighbors),length(n_neighbors)==1,
                          is.numeric(random_state),length(random_state)==1)
  if (n > nrow(dep)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_umap() with 'n' <= ",
                nrow(dep), "."), call. = FALSE)
  }
  if (n_neighbors > ncol(dep)){
    stop(paste0("'n_neighbors' argument is not valid.\n", "Run plot_umap() with 'n_neighbors' <= ",
                ncol(dep), "."), call. = FALSE)
  }
  # for shiny indicate can select the third indicate
  mycolData <- colData(dep)
  mycolData$Condition = mycolData$condition
  mycolData$Replicate = mycolData$replicate
  colData(dep) = mycolData

  columns <- colnames(colData(dep))
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_umap() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ", deparse(substitute(dep)),
                  ".\nValid columns are: '", paste(columns,
                                                   collapse = "', '"), "'."), call. = FALSE)
    }
  }
  var <- apply(assay(dep), 1, sd)
  # n = min(n,nrow(assay(dep)))
  # n_neighbors = min( nrow(t(df)),n_neighbors)

  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)], ]
  umap_res = umap::umap(t(df),n_neighbors= n_neighbors,random_state=random_state)
  umap_df  = umap_res$layout %>% data.frame() %>% rownames_to_column() %>%
    left_join(., data.frame(colData(dep)), by = c(rowname = "ID"))
  colnames(umap_df)[2:3] = c("x","y")

  limit = unlist(umap_df[,2:3])
  p = ggplot(umap_df,aes(x,y)) + labs(title = paste0("UMAP plot total ",n, " variables"))+
    coord_fixed() +
    theme_DEP1()
  if(if_square){
    p <- p +
      scale_x_continuous(limits = range(limit)) +
      scale_y_continuous(limits = range(limit))
  }

  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = umap_df[[indicate[1]]]),
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = umap_df[[indicate[1]]],
                            shape = umap_df[[indicate[2]]]), size = point_size) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = umap_df[[indicate[1]]],
                            shape = umap_df[[indicate[2]]]), size = point_size) +
      facet_wrap(~umap_df[[indicate[3]]]) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    p <- p + labs(x = "", y = "")
    return(p)
  }
  else {
    df <- umap_df %>% dplyr::select(rowname, c("x","y"), match(indicate, colnames(umap_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}


plot_Tsne <- function (dep, indicate = c("condition",
                                         "replicate"), label = FALSE, n = nrow(dep), point_size = 4,
                       label_size = 3, plot = TRUE, if_square = FALSE,
                       perplexity = 2, theseed = NULL)
{
  # if (is.integer(x))
  #   x <- as.numeric(x)
  # if (is.integer(y))
  #   y <- as.numeric(y)
  if (is.integer(n))
    n <- as.numeric(n)
  if (is.integer(point_size))
    point_size <- as.numeric(point_size)
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.numeric(n), length(n) == 1, is.character(indicate),
                          is.logical(label), is.numeric(point_size), length(point_size) ==1,
                          is.numeric(label_size), length(label_size) ==1,
                          is.logical(plot), length(plot) == 1,
                          is.numeric(perplexity),length(perplexity)==1)
  if (n > nrow(dep)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_umap() with 'n' <= ",
                nrow(dep), "."), call. = FALSE)
  }
  # if (n_neighbors > ncol(dep)){
  #   stop(paste0("'n_neighbors' argument is not valid.\n", "Run plot_Tsne() with 'n_neighbors' <= ",
  #               ncol(dep), "."), call. = FALSE)
  # }
  # for shiny indicate can select the third indicate
  mycolData <- colData(dep)
  mycolData$Condition = mycolData$condition
  mycolData$Replicate = mycolData$replicate
  colData(dep) = mycolData

  columns <- colnames(colData(dep))
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_Tsne() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ", deparse(substitute(dep)),
                  ".\nValid columns are: '", paste(columns,
                                                   collapse = "', '"), "'."), call. = FALSE)
    }
  }
  var <- apply(assay(dep), 1, sd)
  # n = min(n,nrow(assay(dep)))
  # n_neighbors = min( nrow(t(df)),n_neighbors)

  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)], ]
  if(is.numeric(theseed)&length(theseed)==1)
    set.seed(theseed)

  Tsne_res = Rtsne::Rtsne(t(df),perplexity = perplexity)
  Tsne_df  = Tsne_res$Y %>% data.frame() %>% `rownames<-`(colnames(df)) %>%rownames_to_column() %>%
    left_join(., data.frame(colData(dep)), by = c(rowname = "ID"))
  colnames(Tsne_df)[2:3] = c("x","y")

  limit = unlist(Tsne_df[,2:3])
  p = ggplot(Tsne_df,aes(x,y)) + labs(title = paste0("T-SNE plot total ",n, " variables"))+
    coord_fixed() +
    theme_DEP1()
  if(if_square){
    p <- p +
      scale_x_continuous(limits = range(limit)) +
      scale_y_continuous(limits = range(limit))
  }

  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = Tsne_df[[indicate[1]]]),
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = Tsne_df[[indicate[1]]],
                            shape = Tsne_df[[indicate[2]]]), size = point_size) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = Tsne_df[[indicate[1]]],
                            shape = Tsne_df[[indicate[2]]]), size = point_size) +
      facet_wrap(~Tsne_df[[indicate[3]]]) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    p <- p + labs(x = "", y = "")
    return(p)
  }
  else {
    df <- Tsne_df %>% dplyr::select(rowname, c("x","y"), match(indicate, colnames(Tsne_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }

}
############


########################################################################### For DEP-RNAseq part, add two plots after pca plot
library(umap)
library(Rtsne)
# /**
#  * [pca plot]
#  * @param  {[type]} data       [row gene , column sample]
#  * @param  {[type]} coldata    [dataframe contain label, ID,condition and replicate,eg:]
#  # DataFrame with 12 rows and 4 columns
# #              label          ID   condition replicate
# #        <character> <character> <character> <numeric>
# # Ubi4_1      Ubi4_1      Ubi4_1        Ubi4         1
# # Ubi4_2      Ubi4_2      Ubi4_2        Ubi4         2
# # Ubi4_3      Ubi4_3      Ubi4_3        Ubi4         3
# # Ubi6_1      Ubi6_1      Ubi6_1        Ubi6         1
# # Ubi6_2      Ubi6_2      Ubi6_2        Ubi6         2
# # ...            ...         ...         ...       ...
# # Ctrl_2      Ctrl_2      Ctrl_2        Ctrl         2
# # Ctrl_3      Ctrl_3      Ctrl_3        Ctrl         3
# # Ubi1_1      Ubi1_1      Ubi1_1        Ubi1         1
# # Ubi1_2      Ubi1_2      Ubi1_2        Ubi1         2
# # Ubi1_3      Ubi1_3      Ubi1_3        Ubi1         3
#  * @param  {[type]} indicate   [Character, Sets the color, shape and facet_wrap of the plot based on columns from the experimental coldata]
#  * @param  {[type]} label      [Logical, Whether or not to add sample labels]
#  * @param  {Number} n          [Integer(1), Sets the number of top variable proteins to consider, default all rows]
#  * @param  {Number} point_size [Integer(1), Sets the size of the points]
#  * @param  {Number} label_size [Integer(1), Sets the size of the labels]
#  * @param  {[type]} plot       [Logical(1), If TRUE (default) the PCA plot is produced. Otherwise (if FALSE), the data which the PCA plot is based on are returned]
#  * @param  {[type]} if_square  [if xlim == ylim, default FALSE]
#  * @param  {[type]} n_neighbors  [integer; number of nearest neighbors]
#  * @param  {[type]} random_state  [ integer; seed for random number generation used during predict()]

#  * @return {[type]}            [description]
#  */
plot_umap_for_rnaseq <- function (data, coldata, indicate = c("condition",
                                                              "replicate"), label = FALSE, n = nrow(data), point_size = 4,
                                  label_size = 3, plot = TRUE, if_square = FALSE,
                                  n_neighbors=4,random_state=10086)
{
  if (is.integer(n))
    n <- as.numeric(n)
  if (is.integer(point_size))
    point_size <- as.numeric(point_size)
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)
  assertthat::assert_that(is.numeric(n), length(n) == 1, is.character(indicate),
                          is.logical(label), is.numeric(point_size), length(point_size) ==1,
                          is.numeric(label_size), length(label_size) ==1,
                          is.logical(plot), length(plot) == 1,
                          is.numeric(n_neighbors),length(n_neighbors)==1,
                          is.numeric(random_state),length(random_state)==1)


  if (n > nrow(data)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_umap_for_rnaseq() with 'n' <= ",
                nrow(data), "."), call. = FALSE)
  }
  if (n_neighbors > ncol(data)){
    stop(paste0("'n_neighbors' argument is not valid.\n", "Run plot_umap_for_rnaseq() with 'n_neighbors' <= ",
                ncol(data), "."), call. = FALSE)
  }
  # for shiny indicate can select the third indicate
  mycolData <- coldata
  mycolData$Condition = mycolData$condition
  mycolData$Replicate = mycolData$replicate
  columns <- colnames(mycolData)

  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_umap_for_rnaseq() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ", deparse(substitute(data)),
                  ".\nValid columns are: '", paste(columns,
                                                   collapse = "', '"), "'."), call. = FALSE)
    }
  }

  var <- apply(data, 1, sd)
  # n = min(n,nrow(assay(dep)))
  # n_neighbors = min( nrow(t(df)),n_neighbors)
  df <- data[order(var, decreasing = TRUE)[seq_len(n)], ]
  umap_res = umap::umap(t(df),n_neighbors= n_neighbors,random_state=random_state)
  umap_df  = umap_res$layout %>% data.frame() %>% rownames_to_column() %>%
    left_join(., data.frame(mycolData), by = c(rowname = "ID"))
  colnames(umap_df)[2:3] = c("x","y")

  limit = unlist(umap_df[,2:3])
  p = ggplot(umap_df,aes(x,y)) + labs(title = paste0("UMAP plot total ",n, " variables"))+
    coord_fixed() +
    theme_DEP1()
  if(if_square){
    p <- p +
      scale_x_continuous(limits = range(limit)) +
      scale_y_continuous(limits = range(limit))
  }

  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = umap_df[[indicate[1]]]),
                        size = point_size) + labs(col = as.factor(indicate[1]))
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = as.factor(umap_df[[indicate[1]]]),
                            shape = as.factor(umap_df[[indicate[2]]])), size = point_size) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = as.factor(umap_df[[indicate[1]]]),
                            shape = as.factor(umap_df[[indicate[2]]])), size = point_size) +
      facet_wrap(~umap_df[[indicate[3]]]) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    p <- p + labs(x = "", y = "")
    return(p)
  }
  else {
    df <- umap_df %>% dplyr::select(rowname, c("x","y"), match(indicate, colnames(umap_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}
# /**
#  * [pca plot]
#  * @param  {[type]} data       [row gene , column sample]
#  * @param  {[type]} coldata    [dataframe contain label, ID,condition and replicate,eg:]
#  # DataFrame with 12 rows and 4 columns
# #              label          ID   condition replicate
# #        <character> <character> <character> <numeric>
# # Ubi4_1      Ubi4_1      Ubi4_1        Ubi4         1
# # Ubi4_2      Ubi4_2      Ubi4_2        Ubi4         2
# # Ubi4_3      Ubi4_3      Ubi4_3        Ubi4         3
# # Ubi6_1      Ubi6_1      Ubi6_1        Ubi6         1
# # Ubi6_2      Ubi6_2      Ubi6_2        Ubi6         2
# # ...            ...         ...         ...       ...
# # Ctrl_2      Ctrl_2      Ctrl_2        Ctrl         2
# # Ctrl_3      Ctrl_3      Ctrl_3        Ctrl         3
# # Ubi1_1      Ubi1_1      Ubi1_1        Ubi1         1
# # Ubi1_2      Ubi1_2      Ubi1_2        Ubi1         2
# # Ubi1_3      Ubi1_3      Ubi1_3        Ubi1         3
#  * @param  {[type]} indicate   [Character, Sets the color, shape and facet_wrap of the plot based on columns from the experimental coldata]
#  * @param  {[type]} label      [Logical, Whether or not to add sample labels]
#  * @param  {Number} n          [Integer(1), Sets the number of top variable proteins to consider, default all rows]
#  * @param  {Number} point_size [Integer(1), Sets the size of the points]
#  * @param  {Number} label_size [Integer(1), Sets the size of the labels]
#  * @param  {[type]} plot       [Logical(1), If TRUE (default) the PCA plot is produced. Otherwise (if FALSE), the data which the PCA plot is based on are returned]
#  * @param  {[type]} if_square  [if xlim == ylim, default FALSE]
#  * @param  {[type]} perplexity  [numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1]
#  * @param  {[type]} theseed  [ set seed for T-SNE]

#  * @return {[type]}            [description]
#  */
plot_Tsne_for_rnaseq <- function (data, coldata, indicate = c("condition",
                                         "replicate"), label = FALSE, n = nrow(data), point_size = 4,
                       label_size = 3, plot = TRUE, if_square = FALSE,
                       perplexity = 2, theseed = NULL)
{
  # if (is.integer(x))
  #   x <- as.numeric(x)
  # if (is.integer(y))
  #   y <- as.numeric(y)
  if (is.integer(n))
    n <- as.numeric(n)
  if (is.integer(point_size))
    point_size <- as.numeric(point_size)
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)
  assertthat::assert_that(is.numeric(n), length(n) == 1, is.character(indicate),
                          is.logical(label), is.numeric(point_size), length(point_size) ==1,
                          is.numeric(label_size), length(label_size) ==1,
                          is.logical(plot), length(plot) == 1,
                          is.numeric(perplexity),length(perplexity)==1)
  if (n > nrow(data)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_Tsne_for_rnaseq() with 'n' <= ",
                nrow(data), "."), call. = FALSE)
  }
  # if (n_neighbors > ncol(dep)){
  #   stop(paste0("'n_neighbors' argument is not valid.\n", "Run plot_Tsne() with 'n_neighbors' <= ",
  #               ncol(dep), "."), call. = FALSE)
  # }
  # for shiny indicate can select the third indicate
  mycolData <- coldata
  mycolData$Condition = mycolData$condition
  mycolData$Replicate = mycolData$replicate
  columns <- colnames(mycolData)

  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_Tsne_for_rnaseq() with a maximum of 3 indicate features")
    }
    if (any(!indicate %in% columns)) {
      stop(paste0("'", paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ", deparse(substitute(data)),
                  ".\nValid columns are: '", paste(columns,
                                                   collapse = "', '"), "'."), call. = FALSE)
    }
  }
  var <- apply(data, 1, sd)
  # n = min(n,nrow(assay(dep)))
  # n_neighbors = min( nrow(t(df)),n_neighbors)
  df <- data[order(var, decreasing = TRUE)[seq_len(n)], ]

  if(is.numeric(theseed)&length(theseed)==1)
    set.seed(theseed)

  Tsne_res = Rtsne::Rtsne(t(df),perplexity = perplexity)
  Tsne_df  = Tsne_res$Y %>% data.frame() %>% `rownames<-`(colnames(df)) %>%rownames_to_column() %>%
    left_join(., data.frame(mycolData), by = c(rowname = "ID"))
  colnames(Tsne_df)[2:3] = c("x","y")

  limit = unlist(Tsne_df[,2:3])
  p = ggplot(Tsne_df,aes(x,y)) + labs(title = paste0("T-SNE plot total ",n, " variables"))+
    coord_fixed() +
    theme_DEP1()
  if(if_square){
    p <- p +
      scale_x_continuous(limits = range(limit)) +
      scale_y_continuous(limits = range(limit))
  }

  if (length(indicate) == 0) {
    p <- p + geom_point(size = point_size)
  }
  if (length(indicate) == 1) {
    p <- p + geom_point(aes(col = as.factor(Tsne_df[[indicate[1]]])),
                        size = point_size) + labs(col = indicate[1])
  }
  if (length(indicate) == 2) {
    p <- p + geom_point(aes(col = as.factor(Tsne_df[[indicate[1]]]),
                            shape = as.factor(Tsne_df[[indicate[2]]])), size = point_size) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (length(indicate) == 3) {
    p <- p + geom_point(aes(col = as.factor(Tsne_df[[indicate[1]]]),
                            shape = as.factor(Tsne_df[[indicate[2]]])), size = point_size) +
      facet_wrap(~Tsne_df[[indicate[3]]]) +
      labs(col = indicate[1], shape = indicate[2])
  }
  if (label) {
    p <- p + geom_text(aes(label = rowname), size = label_size)
  }
  if (plot) {
    p <- p + labs(x = "", y = "")
    return(p)
  }
  else {
    df <- Tsne_df %>% dplyr::select(rowname, c("x","y"), match(indicate, colnames(Tsne_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }

}


############

#' \code{\link{make_se_parse}()}). It is adviced to first remove
#' proteins with too many missing values using \code{\link{filter_missval}()},
#' normalize the data using \code{\link{normalize_vsn}()} and
#' impute remaining missing values using \code{\link{impute}()}.
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
#' These should be formatted as "SampleA_vs_SampleB" or
#' c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @param design_formula Formula,
#' Used to create the design matrix.
#' @return A SummarizedExperiment object
#' containing fdr estimates of differential expression.
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
#' @export
# se = my_imp;type="control";control="wt"; fdr.type = "BH"
test_diff <- function(se, type = c("control", "all", "manual"),
                      control = NULL, test = NULL,
                      design_formula = formula(~ 0 + condition),
                      # correct.p = T,
                      # p.correct.stastic = c("normal.t","student.t","p.val"),
                      # pval.type = c("raw","corrected"),
                      # p.adjusted.method = c("fdrtool","BH"),
                      # cutoff.method = c("fndr","pct0","locfdr"),
                      fdr.type = c("Strimmer's qvalue(t)","Strimmer's qvalue(p)","BH","Storey's qvalue")) {

  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")

  # Show error if inputs do not contain required columns
  type <- match.arg(type)

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
  fdr.type <- match.arg(fdr.type)
  cat(fdr.type)
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
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit, fdr.type){
    res <- topTable(fit, sort.by = "t", coef = comp,
                    number = Inf, confint = TRUE)
    res <- res[!is.na(res$t),]

    if(fdr.type == "Strimmer's qvalue(t)"){
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
    # fdr_res_saved <<- fdr_res
    # res$qval <- fdr_res$qval  ##classic tail area-based FDR
    # res$lfdr <- fdr_res$lfdr  ## local FDR
    # res$pBH <- p.adjust(fdr_res$pval,method = "BH")
    # if(p.correct.stastic!= "p.val" & pval.type == "corrected") res$P.Value <- fdr_res$pval ## replace the p value from limma by the fdrtool's p
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }

  # Retrieve the differential expression test results
  limma_res <- purrr::map_df(cntrst, retrieve_fun,fdr.type = fdr.type)

  table <- limma_res %>%
    select(rowname, logFC, CI.L, CI.R, t,P.Value, qval, comparison) %>%
    mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    gather(variable, value, -c(rowname,comparison)) %>%
    mutate(variable = recode(variable, logFC = "diff", t="t.stastic", P.Value = "p.val", qval = "p.adj")) %>%
    unite(temp, comparison, variable) %>%
    spread(temp, value)

  # # Select the logFC, CI and qval variables
  # if(fdr.type=="qval"){
  #   table <- limma_res %>%
  #     select(rowname, logFC, CI.L, CI.R, t,P.Value, qval, comparison) %>%
  #     mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
  #     gather(variable, value, -c(rowname,comparison)) %>%
  #     mutate(variable = recode(variable, logFC = "diff", t="t.stastic", P.Value = "p.val", qval = "p.adj")) %>%
  #     unite(temp, comparison, variable) %>%
  #     spread(temp, value)
  # }else if(fdr.type=="lfdr"){
  #   table <- limma_res %>%
  #     select(rowname, logFC, CI.L, CI.R, t, P.Value, lfdr, comparison) %>%
  #     mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
  #     gather(variable, value, -c(rowname,comparison)) %>%
  #     mutate(variable = recode(variable, logFC = "diff", t="t.stastic", P.Value = "p.val", lfdr = "p.adj")) %>%
  #     unite(temp, comparison, variable) %>%
  #     spread(temp, value)
  # }else if(fdr.type == "BH"){
  #   table <- limma_res %>%
  #     select(rowname, logFC, CI.L, CI.R, t, P.Value, pBH, comparison) %>%
  #     mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
  #     gather(variable, value, -c(rowname,comparison)) %>%
  #     mutate(variable = recode(variable, logFC = "diff", t="t.stastic", P.Value = "p.val", pBH = "p.adj")) %>%
  #     unite(temp, comparison, variable) %>%
  #     spread(temp, value)
  # }
  table_save <<- table

  # if(p.adjusted.method=="BH"){
  #   table$p.adj = p.adjust(table$p.val,method = "BH")
  # }
  rowData(se) <- merge(rowData(se, use.names = FALSE), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE, sort=FALSE)
  return(se)
}
environment(test_diff) = asNamespace("DEP2")


fdrtool <- function (x, statistic = c("normal", "correlation", "pvalue","studentt"),
                     plot = TRUE, color.figure = TRUE, verbose = TRUE, cutoff.method = c("fndr",
                                                                                         "pct0", "locfdr"), pct0 = 0.75)
{
  statistic = match.arg(statistic)
  cutoff.method = match.arg(cutoff.method)
  if (is.vector(x) == FALSE)
    stop("input test statistics must be given as a vector!")
  if (length(x) < 200)
    warning("There may be too few input test statistics for reliable FDR calculations!")
  if (statistic == "pvalue") {
    if (max(x) > 1 | min(x) < 0)
      stop("input p-values must all be in the range 0 to 1!")
  }
  if (verbose)
    cat("Step 1... determine cutoff point\n")
  if (cutoff.method == "pct0") {
    if (statistic == "pvalue")
      x0 = quantile(x, probs = 1 - pct0)
    else x0 = quantile(abs(x), probs = pct0)
  }
  else if (cutoff.method == "locfdr" & (statistic == "normal" |
                                        statistic == "correlation")) {
    if (statistic == "normal")
      z = x
    if (statistic == "correlation")
      z = atanh(x)
    iqr = as.double(diff(quantile(z, probs = c(0.25, 0.75))))
    sdhat = iqr/(2 * qnorm(0.75))
    N = length(z)
    b = ifelse(N > 5e+05, 1, 4.3 * exp(-0.26 * log(N, 10)))
    z0 = b * sdhat
    if (statistic == "normal")
      x0 = z0
    if (statistic == "correlation")
      x0 = tanh(z0)
  }
  else {
    if (cutoff.method == "locfdr")
      warning("cutoff.method=\"locfdr\" only available for normal and correlation statistic.")
    x0 = fndr.cutoff(x, statistic)
  }
  if (verbose)
    cat("Step 2... estimate parameters of null distribution and eta0\n")
  cf.out <- censored.fit(x = x, cutoff = x0, statistic = statistic)
  if (statistic == "pvalue")
    scale.param = NULL
  else scale.param <- cf.out[1, 5]
  eta0 = cf.out[1, 3]
  if (verbose)
    cat("Step 3... compute p-values and estimate empirical PDF/CDF\n")
  nm = get.nullmodel(statistic)
  pval = nm$get.pval(x, scale.param)
  ee <- ecdf.pval(pval, eta0 = eta0)
  g.pval <- grenander(ee)
  f.pval = approxfun(g.pval$x.knots, g.pval$f.knots, method = "constant",
                     rule = 2)
  f0.pval = function(x) return(ifelse(x > 1 | x < 0, 0, rep(1,
                                                            length(x))))
  F.pval = approxfun(g.pval$x.knots, g.pval$F.knots, method = "linear",
                     yleft = 0, yright = g.pval$F.knots[length(g.pval$F.knots)])
  F0.pval = function(x) return(ifelse(x > 1, 1, ifelse(x <
                                                         0, 0, x)))
  fdr.pval = function(p) {
    p[p == .Machine$double.eps] = 0
    pmin(eta0/f.pval(p), 1)
  }
  Fdr.pval = function(p) pmin(eta0 * p/F.pval(p), 1)
  if (verbose)
    cat("Step 4... compute q-values and local fdr\n")
  qval <- Fdr.pval(pval)
  lfdr <- fdr.pval(pval)
  result = list(pval = pval, qval = qval, lfdr = lfdr, statistic = statistic,
                param = cf.out)
  if (plot) {
    if (verbose)
      cat("Step 5... prepare for plotting\n")
    if (statistic == "pvalue") {
      f0 <- function(zeta) return(nm$f0(zeta, scale.param))
      F0 <- function(zeta) return(nm$F0(zeta, scale.param))
      get.pval <- function(zeta) return(nm$get.pval(1 -
                                                      zeta, scale.param))
      x0 = 1 - x0
    }
    else {
      f0 <- function(zeta) return(2 * nm$f0(zeta, scale.param))
      F0 <- function(zeta) return(2 * nm$F0(zeta, scale.param) -
                                    1)
      get.pval <- function(zeta) return(nm$get.pval(zeta,
                                                    scale.param))
    }
    fdr = function(zeta) fdr.pval(get.pval(zeta))
    Fdr = function(zeta) Fdr.pval(get.pval(zeta))
    F = function(zeta) 1 - eta0 * get.pval(zeta)/Fdr(zeta)
    FA = function(zeta) (F(zeta) - eta0 * F0(zeta))/(1 -
                                                       eta0)
    f = function(zeta) eta0 * (f0(zeta))/fdr(zeta)
    fA = function(zeta) (f(zeta) - eta0 * f0(zeta))/(1 -
                                                       eta0)
    ax = abs(x)
    if (statistic == "pvalue")
      ax = 1 - ax
    xxx = seq(0, max(ax), length.out = 500)
    ll = pvt.plotlabels(statistic, scale.param, eta0)
    par(mfrow = c(3, 1))
    if (color.figure)
      cols = c(2, 4)
    else cols = c(1, 1)
    hist(ax, freq = FALSE, bre = 50, main = ll$main, xlab = ll$xlab,
         cex.main = 1.8)
    lines(xxx, eta0 * f0(xxx), col = cols[1], lwd = 2, lty = 3)
    lines(xxx, (1 - eta0) * fA(xxx), col = cols[2], lwd = 2)
    if (statistic == "pvalue")
      pos1 = "topleft"
    else pos1 = "topright"
    legend(pos1, c("Mixture", "Null Component", "Alternative Component"),
           lwd = c(1, 2, 2), col = c(1, cols), lty = c(1, 3,
                                                       1), bty = "n", cex = 1.5)
    plot(xxx, F(xxx), lwd = 1, type = "l", ylim = c(0, 1),
         main = "Density (first row) and Distribution Function (second row)",
         xlab = ll$xlab, ylab = "CDF", cex.main = 1.5)
    lines(xxx, eta0 * F0(xxx), col = cols[1], lwd = 2, lty = 3)
    lines(xxx, (1 - eta0) * FA(xxx), col = cols[2], lwd = 2)
    plot(xxx, Fdr(xxx), type = "l", lwd = 2, ylim = c(0,
                                                      1), main = "(Local) False Discovery Rate", ylab = "Fdr and fdr",
         xlab = ll$xlab, lty = 3, cex.main = 1.5)
    lines(xxx, fdr(xxx), lwd = 2)
    if (eta0 > 0.98)
      pos2 = "bottomleft"
    else pos2 = "topright"
    legend(pos2, c("fdr (density-based)", "Fdr (tail area-based)"),
           lwd = c(2, 2), lty = c(1, 3), bty = "n", cex = 1.5)
    par(mfrow = c(1, 1))
    rm(ax)
  }
  if (verbose)
    cat("\n")
  return(result)
}

environment(fun  = .GlobalEnv$fdrtool) = asNamespace("fdrtool")

# x_saved <<- input$stastical_plot_x; y_saved <<- input$stastical_plot_y
# plottype_saved <<- input$stastical_plot_type
# Volcano_cntrst_stastical_plot_saved <<- input$Volcano_cntrst_stastical_plot
# dep_saved <<- selected()
# dep = dep_saved;x_saved -> x ; y_saved -> y; plottype_saved -> plottype; Volcano_cntrst_stastical_plot_saved -> contrast
statistics_plot <- function(dep,x = c("-log10.padj", "-log10.pval","padj","pval","t") ,y = c("-log10.padj", "-log10.pval","padj","pval","t"), contrast, plottype = c("x-y","histogram") ){
  row_data <- rowData(dep, use.names = FALSE)
  # row_data_saved <<- row_data
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



## plot the fit normal for log2 fold change for DEP-LFQ panel
plot_norm_distribution <- function(dep, contrast) {
  row_data <- as.data.frame(rowData(dep))
  cols_diff <- grep(paste(contrast, "_diff", sep = ""), colnames(row_data))


    x = row_data[, cols_diff]
    x = fun.outlier(x) %>% na.omit()
    fit <- fitnormal(x)
    σ = sqrt(fit$theta[2])
    mu = fit$theta[1]

    df <- as.data.frame(x)
    colnames(df) = gsub("_diff", "", colnames(row_data)[cols_diff])
    dat <- df %>% tibble::rownames_to_column() %>% tidyr::gather("contrast", "diff", -rowname)

    p <- ggplot(dat, aes(x = diff)) +
      geom_histogram(aes(y =..density..),
                     bins = 100,
                     colour ="black",
                     fill ="grey") +
      stat_function(fun = dnorm, args = list(mean = mu, sd = σ), color ="darkred") +
      # annotate("text", x = Inf, y = Inf, label = "Upper left", hjust = 2,
      # vjust = 2) +
      geom_text(data = data.frame(), aes(x = Inf, y = Inf,label = paste("Mean = ", round(mu, 2), "\nStd.Dev. =", round(σ, 2), sep = "")), hjust = 1.2, vjust = 1.2, size = 4) +
      labs(x = "log2 fold change") +
      theme_DEP2() +
      theme(legend.position = "none") +
      facet_wrap(~ contrast)
    return(p)
}

## save functions
#for DEP-LFQ panel,heatmap when if_chooseToshow is TRUE, the choices of proteins
heatmap_name <- function (dep, manual = FALSE, manual_name = NULL, same_trend = FALSE, ...)
{

  # library(tidyr)
  # library(tibble)
  row_data <- rowData(dep, use.names = FALSE)
  # col_data <- colData(dep) %>% as.data.frame()

  # library(stringr)
  if(manual) {
    my_ind <- row_data[ , paste(manual_name, "_significant", sep = "")]
    my = my_ind
    if(class(my_ind) == "DFrame"){
      if(same_trend){
        my1 = apply(my_ind, 1, any)
        #my_ind <- row_data[ , c(paste(manual_name, "_significant", sep = ""), "significant")]
        my_ind = apply(my_ind, 1, function(i){
          #all(i[1:((ncol(my_ind)-1)/2)] == i[((ncol(my_ind)-1 )/2 + 1):(ncol(my_ind)-1)]) & i[ncol(my_ind)]
          all(i[1:((ncol(my_ind))/2)] == i[((ncol(my_ind) )/2 + 1):(ncol(my_ind))])
        })
        ratio <- row_data[ , paste(manual_name, "_diff", sep = "")]
        ratio <- as.matrix(ratio)
        ratio[which(ratio >= 0)] = 1
        ratio[which(ratio < 0)] = -1

        ratio <- apply(ratio, 1, function(myratio){
          all(myratio[1:(ncol(ratio)/2)] == myratio[(ncol(ratio)/2 + 1):ncol(ratio)])
        })

        my_ind = intersect(which(my_ind),which(ratio))
        my_ind  = intersect(my_ind, which(my1))

      } else {
        my_ind = apply(my_ind, 1, any)
      }
    }
    filtered <- dep[my_ind , ]
    #get index of manual_name
    colname = colnames(assay(filtered))
    ind = which(gsub("_[[:digit:]]*", "", colname) %in% unlist(strsplit(manual_name, split = "_vs_")))
  } else {
    filtered <- dep[row_data$significant, ]
  }

  return(rowData(filtered)$name)
}


#for DEG-RNAseq panel,heatmap, return the choices of genes
# heatmap_name_rnaseq <- function (data, row_data, manual = FALSE, manual_name = NULL, ...)
# {
#   row_data <- row_data
#
#   library(stringr)
#   if(manual) {
#     my_ind <- row_data[ , paste(manual_name, "_significant", sep = "")]
#     my = my_ind
#     if(class(my_ind) == "data.frame"){
#       my_ind = apply(my_ind, 1, any)
#
#     }
#     filtered <- data[my_ind , ]
#
#     #get index of manual_name
#     colname = colnames(filtered)
#
#     ind = which(gsub("_[[:digit:]]*", "", colname) %in% unlist(strsplit(manual_name, split = "_vs_")))
#   } else {
#     filtered <- data[row_data$significant, ]
#   }
#   data = filtered
#   return(rownames(data))
# }
#
# heatmap_name <- function (data, manual_name = NULL, ...)
# {
#   DEP2::get_signicant(data, contrasts = manual_name)
#   return(rownames(data))
# }

#for pro-rna heatmap together panel, for upload pro and rna data, heatmap when if_chooseToshow_for_pro_rna is TRUE, the choices of genes
pro_rna_heatmap_name <- function(rna, pro, gene, Type = "all") {
  library(dplyr)
  ##### in order to know which gene no rna or pro data, and set gene color to grey which have no data
  pro$name_pro = pro$name
  rna$name_rna = rna$name

  if(Type == "all genes") {
    gene = gene
  }

  if(Type == "both") {
    gene = data.frame(name = intersect(rna$name[na.omit(match(gene$name, rna$name))], pro$name[na.omit(match(gene$name, pro$name))]), stringsAsFactors = F)
  }

  if(Type == "at least one") {
    gene = data.frame(name = union(rna$name[na.omit(match(gene$name, rna$name))], pro$name[na.omit(match(gene$name, pro$name))]), stringsAsFactors = F)
  }

  df_rna = rna[match(gene$name, rna$name),]
  df_rna$gene = gene$name
  df <- merge(df_rna, pro, by.x = "gene", by.y = "name", all.x = TRUE)
  #rna, pro set rowname color, no data to grey
  df$row_color_RNA = "black"
  df$row_color_RNA[is.na(df$name_rna)] = "grey"

  df$row_color_pro = "black"
  df$row_color_pro[is.na(df$name_pro)] = "grey"
  colnames(df) = gsub("\\.y$", "", gsub("\\.x$", "", colnames(df)))

  dat = df[ , c(match("gene", colnames(df)), match("name_rna", colnames(df)), match("row_color_RNA", colnames(df)), 3:(match("name_rna", colnames(df)) -1), match("name_pro", colnames(df)), match("row_color_pro", colnames(df)), (match("name_rna", colnames(df)) +1):(match("name_pro", colnames(df)) -1))]
  names(dat) = gsub("\\.1","",colnames(dat))

  return(dat$gene)
}



#for pro-rna heatmap together panel, for import from DEP-LFQ and DEG-RNAseq, heatmap when if_chooseToshow_for_pro_rna is TRUE, the choices of genes
pro_rna_heatmap_name_for_whole <- function(rna, pro, manual = FALSE, contrast = NULL, Type = "at least one") {
  library(dplyr)
  library(circlize)
  library(RColorBrewer)
  library(ComplexHeatmap)

  if(manual) {
  gene_rna = rna %>% dplyr::filter(get(paste("rna_", contrast, "_significant", sep = ""))) %>% dplyr::select(rna_name)
  gene_pro = pro %>% dplyr::filter(get(paste("pro_", contrast, "_significant", sep = ""))) %>% dplyr::select(pro_name)
  if(Type == "at least one"){
    gene = data.frame(name = na.omit(union(gene_rna$rna_name, gene_pro$pro_name)))
  } else {
    if(Type == "both significant") {
      gene = data.frame(name = na.omit(intersect(gene_rna$rna_name, gene_pro$pro_name)))
      } else {
        if(Type == "both") {
          gene = na.omit(union(gene_rna$rna_name, gene_pro$pro_name))
          gene = data.frame(name = intersect(rna$rna_name[na.omit(match(gene, rna$rna_name))], pro$pro_name[na.omit(match(gene, pro$pro_name))]), stringsAsFactors = F)
        }
      }
  }
}


  ###### any significant
  if(!manual) {
    gene_rna = rna %>% dplyr::filter(rna_significant) %>% dplyr::select(rna_name)
    gene_pro = pro %>% dplyr::filter(pro_significant) %>% dplyr::select(pro_name)
    if(Type == "at least one") {
      gene = data.frame(name = na.omit(union(gene_rna$rna_name, gene_pro$pro_name)))
    } else {
      if(Type == "both significant") {
        gene = data.frame(name = na.omit(intersect(gene_rna$rna_name, gene_pro$pro_name)))
        } else {
          if(Type == "both") {
            gene = na.omit(union(gene_rna$rna_name, gene_pro$pro_name))
            gene = data.frame(name = intersect(rna$rna_name[na.omit(match(gene, rna$rna_name))], pro$pro_name[na.omit(match(gene, pro$pro_name))]), stringsAsFactors = F)
          }
        }
    }
  }

  ##### in order to know which gene no rna or pro data, and set gene color to grey which have no data
  pro$name_pro = pro$pro_name
  rna$name_rna = rna$rna_name

  df_rna = rna[match(gene$name, rna$rna_name),]
  df_rna$gene = gene$name
  df <- merge(df_rna, pro, by.x = "gene", by.y = "pro_name", all.x = TRUE)

  ###### for contrast
  #rna, pro set rowname color, no data to grey
  if(manual) {
    df$row_color_RNA = "black"
    df$row_color_RNA[df[ , paste("rna_", contrast, "_significant", sep = "")]] = "red"
    df$row_color_RNA[is.na(df$name_rna)] = "grey"

    df$row_color_pro = "black"
    df$row_color_pro[df[ , paste("pro_", contrast, "_significant", sep = "")]] = "red"
    df$row_color_pro[is.na(df$name_pro)] = "grey"

    colnames(df) = gsub("rna_", "", gsub("pro_", "", colnames(df)))

    df = df[ , -grep("_vs_|significant", names(df))]
    names(df) = gsub("\\.1","",colnames(df))
    df = df[ , which(gsub("_\\d*$", "", colnames(df)) %in% c("gene", "name_rna", "name_pro", "row_color_RNA", "row_color_pro", unlist(strsplit(contrast, split = "_vs_"))))]
    names(df) = gsub("\\.1","",colnames(df))
    dat = df[ , c(match("gene", colnames(df)), match("name_rna", colnames(df)), match("row_color_RNA", colnames(df)), 2:(match("name_rna", colnames(df)) -1), match("name_pro", colnames(df)), match("row_color_pro", colnames(df)), (match("name_rna", colnames(df)) +1):(match("name_pro", colnames(df)) -1))]
    names(dat) = gsub("\\.1","",colnames(dat))
  }



  ##### for any significant
  #rna, pro set rowname color, no data to grey
  if(!manual) {
    df$row_color_RNA = "black"
    df$row_color_RNA[df[ , "rna_significant"]] = "red"
    df$row_color_RNA[is.na(df$name_rna)] = "grey"

    df$row_color_pro = "black"
    df$row_color_pro[df[ , "pro_significant"]] = "red"
    df$row_color_pro[is.na(df$name_pro)] = "grey"

    colnames(df) = gsub("rna_", "", gsub("pro_", "", colnames(df)))

    df = df[ , -grep("_vs_|significant", names(df))]
    dat = df[ , c(match("gene", colnames(df)), match("name_rna", colnames(df)), match("row_color_RNA", colnames(df)), 3:(match("name_rna", colnames(df)) -1), match("name_pro", colnames(df)), match("row_color_pro", colnames(df)), (match("name_rna", colnames(df)) +1):(match("name_pro", colnames(df)) -1))]
    names(dat) = gsub("\\.1","",colnames(dat))
  }

  return(dat$gene)
  }

#$add
#dep:  SummarizedExperiment, Data object for which differentially enriched proteins are annotated (output from test_diff() and add_rejections()).
#contrast: Character(1), Specifies the contrast to plot.
#adjusted: Logical(1), Whether or not to use adjusted p values
#labelWay: the way of adding label, one of c("all significant", "up","down","selected proteins"), default: "all significant"
#showNum: the number of the added labels corresponding to the labelWay, default: 20, when showNum = 0, there is no label added, and when it is more than the labelWay, it is the nrow(labelWay), eg: you can set it 10000 in order to show all label of the labelWay
#chooseTolabel: the protein names of your choose , when labelWay == "selected proteins" act , default: NULL
#fontSize: the font size of the label, default: 7
#dotsize: the dot size of the label, default: 0.8
#fcCutoff: log2 fold change cutoff
#adjpCutoff: adj pvalue cutoff
#plot: Logical(1), If TRUE (default) the volcano plot is produced. Otherwise (if FALSE), the data which the volcano plot is based on are returned.
#label.rectangle: TRUE or FALSE, whether rectangle around label
#stroke: label point lwd, default : 0.6
#PTM: if for PTM volcano, default FALSE
my_plot_volcano <- function(dep, contrast, adjusted = FALSE ,labelWay = "all significant",
                            showNum = 20, chooseTolabel = NULL,fontSize = 9, dotsize = 2, same_width = TRUE, plot = TRUE,
                            fcCutoff = 1, linetype = 4, adjpCutoff = 0.05,
                            if_peptide_color = FALSE, label.rectangle = FALSE, stroke = 0.6,
                            down_color = "#0099CC", stable_color = "#BEBEBE", up_color = "#CD8500", peptide_1_color = "#A020F0", peptide_2_color = "#0000FF",
                            threshold_method = c("intersect","curve"),curvature=0.6,x0_fold = 2, PTM = FALSE) {
  library(ggpubr)
  library(ggthemes)
  #  if (is.integer(label_size))
  #     label_size <- as.numeric(label_size)
  # assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
  #     is.character(contrast), length(contrast) == 1, is.numeric(label_size),
  #     length(label_size) == 1, is.logical(add_names), length(add_names) ==
  #         1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot),
  #     length(plot) == 1)
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."),
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) <
      1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if (length(grep("_significant", colnames(row_data))) <
      1) {
    stop(paste0("'[contrast]_significant' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."),
         call. = FALSE)
  }
  if (length(grep(paste(contrast, "_diff", sep = ""),
                  colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>% data.frame() %>% dplyr::select(ends_with("_diff")) %>%
      colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg, call. = FALSE)
  }

  if( length(threshold_method)!=1 | (!threshold_method %in% c("intersect","curve")) ){
    stop("threshold_method should be one of 'intersect'/'curve'",call. = F)
  }
  #when name have Disadvantages eg: when contrast = c("H2A_vs_Biotin","uH2A_vs_Biotin"), when set contrast = c("H2A_vs_Biotin"), can find two cols,and then error, the same problem with p_values and signif
  #diff <- grep(paste(contrast, "_diff", sep = ""),
  #    colnames(row_data))
  diff <- match(paste(contrast, "_diff", sep = ""),
                colnames(row_data))
  if (adjusted) {
    p_values <- match(paste(contrast, "_p.adj", sep = ""),
                      colnames(row_data))
  } else {
    p_values <- match(paste(contrast, "_p.val", sep = ""),
                      colnames(row_data))
  }
  signif <- match(paste(contrast, "_significant", sep = ""),
                  colnames(row_data))

  if(PTM) {
    df <- data.frame(x = row_data[, diff], y = -log10(row_data[, p_values]), significant = row_data[, signif], name = row_data$name) %>% dplyr::filter(!is.na(significant)) %>% arrange(significant)
  } else {
    df <- data.frame(x = row_data[, diff], y = -log10(row_data[, p_values]), significant = row_data[, signif], name = row_data$name, Peptides = row_data$Peptides) %>% dplyr::filter(!is.na(significant)) %>% arrange(significant)
  }

  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  xlimit <- ceiling(max(c(abs(min(df$x)), abs(max(df$x)))))
  ylimit <- ceiling(abs(max(df$y)))


  test <- df
  test$change = ifelse((test$significant),
                       ifelse(test$x > 0 ,'up','down'),
                       'stable')
  #add Peptides of significant proteins, 1 to purple, 2 to blue
  if(if_peptide_color){
    test[!test$change == "stable",]$change[which(test[!test$change == "stable",]$Peptides == 1)] = "peptide_1"
    test[!test$change == "stable",]$change[which(test[!test$change == "stable",]$Peptides == 2)] = "peptide_2"
  }
  # test <- test[order(test$pval,decreasing = F),]
  if(labelWay == "selected proteins"){
    labels <- chooseTolabel
    for_label <<- test[which(test$name %in% labels),]
    if(is.null(labels) | nrow(for_label) == 0) return(print("please choose at least one protein"))
    ########################################################到这里
  }else if(labelWay == "all significant"){
    labelProtein <- test[which(test$change!='stable'),]
    labelProtein <- labelProtein[order(labelProtein$y,decreasing = T),]
    labels <- labelProtein$name[c(0:min(nrow(labelProtein),showNum))]
    for_label <- test[which(test$name %in% labels),]
  }else{
    labelProtein <-test[which(test$change == labelWay),]
    labelProtein <- labelProtein[order(labelProtein$y,decreasing = T),]
    labels <- labelProtein$name[c(0:min(nrow(labelProtein),showNum))]
    for_label <- test[which(test$name %in% labels),]
  }

  if(adjusted) {
    ylab = expression(-log[10] ~ "Adjusted p-value")
  } else {
    ylab = expression(-log[10] ~ "P-value")
  }


  volcanoPlot <-  ggscatter(test,
                            x = "x",
                            y = "y",
                            title = contrast,
                            xlab = expression(log[2] ~ "Fold change"),
                            ylab = ylab,
                            shape = 16,#20
                            color = "change",
                            size = dotsize,
                            label = "name",
                            font.label = c(fontSize, "plain") ,
                            repel = T,
                            palette = c("down" = down_color, 'stable' = stable_color, "up" = up_color, "peptide_1" = ifelse(if_peptide_color, peptide_1_color, "#A020F0"), "peptide_2" = ifelse(if_peptide_color, peptide_2_color, "#0000FF")) ,# #CC9933,#FC4E07,lightcoral,#0099CC,#00AFBB
                            label.select = labels,
                            label.rectangle = label.rectangle) +
    theme_base() +
    #adjust pvalue volcano in DEP is very ugly, Basically not used, so no need to add the limit line based on padj
    # geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col="grey",lwd=0.6) +
    ##### geom_vline(xintercept=c(-fcCutoff,fcCutoff),linetype=linetype, col=c(down_color, up_color),lwd=0.6) + ### 正交的
    theme(plot.title=element_text(hjust=0.5)) +
    theme_DEP1() +
    theme(legend.position = "none") +
    # geom_vline(xintercept = 0) +
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1,
                                                                                      0), vjust = c(-1, -1), label = c(name1, name2), size = 5, fontface = "bold"))+
    # scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, ylimit)) +#limits = c(0, ylimit) expand = expand_scale(mult = 0, add = c(0, 0.15))
    geom_point(aes(fill = change), data = for_label, size = if(labelWay == "selected proteins"){dotsize + 1} else {dotsize}, shape = 21, color = "black", stroke = stroke, position = "identity")

  if(same_width){
    volcanoPlot <- volcanoPlot + scale_x_continuous(limits = c(-xlimit, xlimit))
  }

  if(threshold_method == "intersect"){
    volcanoPlot <- volcanoPlot + geom_vline(xintercept=c(-fcCutoff,fcCutoff),linetype=linetype, col=c(down_color, up_color),lwd=0.6)
    if(adjusted){#for adjust pvalue
      volcanoPlot <- volcanoPlot + geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col=stable_color,lwd=0.6)
    }else {#for pvalue
      volcanoPlot <- volcanoPlot + geom_hline(yintercept=min(test$y[!test$change == "stable"]),linetype=linetype, col=stable_color,lwd=0.6)
      #     volcanoPlot <- volcanoPlot + geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col="grey",lwd=0.6)
    }
  }else if(threshold_method == "curve"){

    x = df[,"x"]
    x = fun.outlier(x) %>% na.omit()
    fit <- fitnormal(x)
    σ = sqrt(fit$theta[2])
    x0 = x0_fold*σ

    curve_fun1 = function(x){
      y=abs(curvature/(x-x0))
    }
    curve_fun2 = function(x){
      y= abs(curvature/(x+x0))
    }
    volcanoPlot <- volcanoPlot + stat_function(fun=curve_fun1,
                                               xlim=c(x0,xlimit),
                                               colour=up_color,linetype=linetype,lwd=0.6)
    volcanoPlot <- volcanoPlot + stat_function(fun=curve_fun2,
                                               xlim=c(-xlimit,-x0),
                                               colour=down_color,linetype=linetype,lwd=0.6)

  }


  if (plot) {
    return(volcanoPlot)
  }
  else {
    test <- test %>% dplyr::select(name, x, y, significant, change) %>% arrange(desc(x))
    colnames(test)[c(1, 2, 3)] <- c("protein", "log2_fold_change",
                                    "p_value_-log10")
    if (adjusted) {
      colnames(test)[3] <- "adjusted_p_value_-log10"
    }
    return(test)
  }
}

#$add
#indicate can be one , two, three element, the third means facet_wrap
# if_square: if xlim == ylim, default FALSE
# name: one of c("proteins" , "peptides"), only for the title of the plot
plot_pca <- function (dep, x = 1, y = 2, indicate = c("condition",
    "replicate"), label = FALSE, n = 500, point_size = 4,
    label_size = 3, plot = TRUE, if_square = FALSE, name = "proteins")
{
    if (is.integer(x))
        x <- as.numeric(x)
    if (is.integer(y))
        y <- as.numeric(y)
    if (is.integer(n))
        n <- as.numeric(n)
    if (is.integer(point_size))
        point_size <- as.numeric(point_size)
    if (is.integer(label_size))
        label_size <- as.numeric(label_size)
    assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
        is.numeric(x), length(x) == 1, is.numeric(y), length(y) ==
            1, is.numeric(n), length(n) == 1, is.character(indicate),
        is.logical(label), is.numeric(point_size), length(point_size) ==
            1, is.numeric(label_size), length(label_size) ==
            1, is.logical(plot), length(plot) == 1)
    if (x > ncol(dep) | y > ncol(dep)) {
        stop(paste0("'x' and/or 'y' arguments are not valid\n",
            "Run plot_pca() with 'x' and 'y' <= ", ncol(dep),
            "."), call. = FALSE)
    }
    if (n > nrow(dep)) {
        stop(paste0("'n' argument is not valid.\n", "Run plot_pca() with 'n' <= ",
            nrow(dep), "."), call. = FALSE)
    }
    # for shiny indicate can select the third indicate
    mycolData <- colData(dep)
    mycolData$Condition = mycolData$condition
    mycolData$Replicate = mycolData$replicate
    colData(dep) = mycolData

    columns <- colnames(colData(dep))
    if (!is.null(indicate)) {
        if (length(indicate) > 3) {
            stop("Too many features in 'indicate'\n        Run plot_pca() with a maximum of 3 indicate features")
        }
        if (any(!indicate %in% columns)) {
            stop(paste0("'", paste0(indicate, collapse = "' and/or '"),
                "' column(s) is/are not present in ", deparse(substitute(dep)),
                ".\nValid columns are: '", paste(columns,
                  collapse = "', '"), "'."), call. = FALSE)
        }
    }
    var <- apply(assay(dep), 1, sd)
    df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)],
        ]
    pca <- prcomp(t(df), scale = FALSE)
    pca_df <- pca$x %>% data.frame() %>% rownames_to_column() %>%
        left_join(., data.frame(colData(dep)), by = c(rowname = "ID"))
    percent <- round(100 * pca$sdev^2/sum(pca$sdev^2), 1)
    for (feat in indicate) {
        pca_df[[feat]] <- as.factor(pca_df[[feat]])
    }
    limit = unlist(pca_df[,c(paste0("PC", x), paste0("PC", y))])
    p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC",
        y)))) + labs(title = paste0("PCA plot - top ",
       n, " variable ", name), x = paste0("PC",
        x, ": ", percent[x], "%"), y = paste0("PC",
        y, ": ", percent[y], "%")) +
        coord_fixed() +
        theme_DEP1()
    if(if_square){
      p <- p +
        scale_x_continuous(limits = range(limit)) +
        scale_y_continuous(limits = range(limit))
    }


    if (length(indicate) == 0) {
        p <- p + geom_point(size = point_size)
    }
    if (length(indicate) == 1) {
        p <- p + geom_point(aes(col = pca_df[[indicate[1]]]),
            size = point_size) + labs(col = indicate[1])
    }
    if (length(indicate) == 2) {
        p <- p + geom_point(aes(col = pca_df[[indicate[1]]],
            shape = pca_df[[indicate[2]]]), size = point_size) +
            labs(col = indicate[1], shape = indicate[2])
    }
    if (length(indicate) == 3) {
        p <- p + geom_point(aes(col = pca_df[[indicate[1]]],
            shape = pca_df[[indicate[2]]]), size = point_size) +
            facet_wrap(~pca_df[[indicate[3]]]) +
        labs(col = indicate[1], shape = indicate[2])
    }
    if (label) {
        p <- p + geom_text(aes(label = rowname), size = label_size)
    }
    if (plot) {
        return(p)
    }
    else {
        df <- pca_df %>% dplyr::select(rowname, paste0("PC", c(x,
            y)), match(indicate, colnames(pca_df)))
        colnames(df)[1] <- "sample"
        return(df)
    }
}



#$add
## name: one of c("proteins" , "peptides"), only for the title of the plot , default: proteins
plot_numbers <- function (se, plot = TRUE, name = "proteins")
{
    assertthat::assert_that(inherits(se, "SummarizedExperiment"),
        is.logical(plot), length(plot) == 1)
    df <- assay(se) %>% data.frame() %>% rownames_to_column() %>%
        gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin),
        0, 1))
    stat <- df %>% group_by(ID) %>% summarize(n = n(), sum = sum(bin)) %>%
        left_join(., data.frame(colData(se)), by = "ID")

    firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
    }

    p <- ggplot(stat, aes(x = ID, y = sum, fill = condition)) +
        geom_col() + geom_hline(yintercept = unique(stat$n)) +
        labs(title = firstup(paste0(name,  " per sample")), x = "",
            y = paste0("Number of ", name)) + theme_DEP2()
    if (plot) {
        return(p)
    }
    else {
        df <- as.data.frame(stat)
        colnames(df)[seq_len(3)] <- c("sample", paste0("total_", name),
            paste0(name, "_in_sample"))
        return(df)
    }
}


#$add
## name: one of c("proteins" , "peptides"), only for the title of the plot , default: proteins
plot_coverage <- function (se, plot = TRUE, name = "proteins")
{
    assertthat::assert_that(inherits(se, "SummarizedExperiment"),
        is.logical(plot), length(plot) == 1)
    df <- assay(se) %>% data.frame() %>% rownames_to_column() %>%
        gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin),
        0, 1))
    stat <- df %>% group_by(rowname) %>% summarize(sum = sum(bin))
    table <- table(stat$sum) %>% data.frame()
    firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
    }

    p <- ggplot(table, aes(x = "all", y = Freq, fill = Var1)) +
        geom_col(col = "white") + scale_fill_grey(start = 0.8,
        end = 0.2) + labs(title = firstup(paste0(substr(name, 1, (nchar(name) -1)), " coverage")), x = "",
        y = paste0("Number of ", name), fill = "Samples") +
        theme_DEP1()
    if (plot) {
        return(p)
    }
    else {
        df <- as.data.frame(table)
        colnames(df) <- c("samples", name)
        return(df)
    }
}



make_se_parse = function (proteins_unique, columns, mode = c("char", "delim"),
                          chars = 1, sep = "_", remove_prefix = T)
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.integer(columns),
                          is.character(mode), is.numeric(chars), length(chars) ==
                            1, is.character(sep), length(sep) == 1)
  mode <- match.arg(mode)
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  if(remove_prefix){
    colnames(raw) <- delete_prefix(colnames(raw)) %>% make.names()
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
    colnames(raw) = gsub(get_suffix(colnames(raw)),"", colnames(raw))
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
environment(make_se_parse) = asNamespace("DEP2")

make_se <- function (proteins_unique, columns, expdesign)
{
  assertthat::assert_that(is.data.frame(proteins_unique), is.integer(columns),
                          is.data.frame(expdesign))
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
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if (tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
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
environment(make_se) = asNamespace("DEP2")

# impute <- function (se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet",
#                                 "MinProb", "man", "min", "zero", "mixed", "nbavg","RF"), ...)
# {
#   assertthat::assert_that(inherits(se, "SummarizedExperiment"),
#                           is.character(fun))
#   fun <- match.arg(fun)
#   if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
#     stop("'name' and/or 'ID' columns are not present in '",
#          deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns",
#          call. = FALSE)
#   }
#   if (!any(is.na(assay(se)))) {
#     warning("No missing values in '", deparse(substitute(se)),
#             "'. ", "Returning the unchanged object.", call. = FALSE)
#     return(se)
#   }
#   rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
#   rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
#   if (fun == "man") {
#     se <- manual_impute(se, ...)
#   }else if(fun == "RF"){
#     doParallel::registerDoParallel(cores=4)
#     assay_save <<- assay(se)
#     assay(se) <- missForest::missForest(assay(se),ntree = 60,parallelize="variables")$ximp
#     library(microbenchmark)
#     # microbenchmark({
#     #   missForest::missForest(assay_save,parallelize="variables")
#     # },times = 8)
#     # microbenchmark({
#     #   missForest::missForest(assay_save,parallelize="forests")
#     # },times = 8)
#     # microbenchmark({
#     #   missForest::missForest(assay_save,ntree = 50,parallelize="variables")
#     # },times = 8)
#
#   }
#   else {
#     MSnSet_data <- as(se, "MSnSet")
#     MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun,
#                                       ...)
#     assay(se) <- MSnbase::exprs(MSnSet_imputed)
#   }
#   return(se)
# }
# environment(fun  = .GlobalEnv$impute) = asNamespace("DEP2")


# ##### $add  for directly upload peptide.txt, and use it assemble to proteingroup.txt, and then do the differential analysis
# ##' normalize_method should be one of MsCoreUtils::normalizeMethods(), detail see QFeatures::normalize
# aggregation <- function(Peptide,
#                         ecols,
#                         fnames,
#                         expdesign = NULL,
#                         NAnum = NULL,
#                         thr = NULL,
#                         normalize_method = c("quantiles","quantiles.robust","vsn","center.mean",
#                                              "center.median","div.mean","div.median", "diff.median"),
#                         aggregate_fun = MsCoreUtils::robustSummary){
#   if (any(!apply(Peptide[, ecols], 2, is.numeric))) {
#     stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
#          call. = FALSE)
#   }
#
#   if (any(!c("label", "condition", "replicate") %in%
#           colnames(expdesign))) {
#     stop("'label', 'condition' and/or 'replicate' columns",
#          "are not present in the experimental design",
#          call. = FALSE)
#   }
#   if (tibble::is_tibble(Peptide))
#     Peptide <- as.data.frame(Peptide)
#   if (tibble::is_tibble(expdesign))
#     expdesign <- as.data.frame(expdesign)
#
#   rownames(expdesign) = expdesign$label
#   Peptide[,ecols] = apply(Peptide[,ecols], 2, function(x){
#     # temp1 <<- x
#     # temp <<- as.integer(x)
#     x[!(!grepl("[A-z]",x) & grepl("\\d",x))] = 0
#     return(as.numeric(x))
#   })
#
#   QF <- readQFeatures(table = Peptide, ecol = ecols,
#                       fnames = fnames,
#                       name = "peptideRaw")
#   rowData(QF[["peptideRaw"]])$nNonZero <- rowSums(assay(QF[["peptideRaw"]]) > 0)
#   # if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))
#
#   QF[["peptideRaw"]]@colData = as(expdesign, "DataFrame")
#   colData(QF) = as(expdesign, "DataFrame")
#   QF <- zeroIsNA(QF, "peptideRaw")
#
#   ## 注：同DEP的DEP::filter_missval
#   filter_thr = function(se,thr){
#     if (any(!c("label", "condition", "replicate") %in%
#             colnames(colData(se)))) {
#       stop("'label', 'condition' and/or 'replicate' columns are not present in '",
#            deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
#            call. = FALSE)
#     }
#     max_repl <- max(colData(se)$replicate)
#     if (thr < 0 | thr > max_repl) {
#       stop("invalid filter threshold applied", "\nRun filter_missval() with a threshold ranging from 0 to ",
#            max_repl)
#     }
#     bin_data <- assay(se)
#     idx <- (is.na(assay(se)))
#     bin_data[!idx] <- 1
#     bin_data[idx] <- 0
#     keep <- bin_data %>% data.frame() %>% rownames_to_column() %>%
#       gather(ID, value, -rowname) %>% left_join(., data.frame(colData(se)),
#                                                 by = "ID") %>% group_by(rowname, condition) %>%
#       summarize(miss_val = n() - sum(value)) %>% filter(miss_val <=
#                                                           thr) %>% spread(condition, miss_val)
#     se_fltrd <- se[keep$rowname, ]
#     return(se_fltrd)
#   }
#
#   if(!is.null(thr))
#     QF[["peptideRaw"]] = filter_thr(QF[["peptideRaw"]], thr = thr)
#
#   ##
#   QF <- logTransform(QF, base = 2, i = "peptideRaw", name = "peptideLog")
#   QF <- QFeatures::normalize(QF,
#                              i = "peptideLog",
#                              name = "peptideNorm",
#                              method = normalize_method)
#   # QF[["peptideNorm"]]
#   limma::plotDensities(assay(QF[["peptideLog"]]))
#   limma::plotDensities(assay(QF[["peptideNorm"]]))#到这里
#   QF <- filterFeatures(QF, ~ Proteins %in% smallestUniqueGroups(rowData(QF[["peptideLog"]])$Proteins))
#   QF_aggregated <- aggregateFeatures(QF,
#                                      i = "peptideNorm", fcol = "Proteins",
#                                      name = "protein",
#                                      fun = aggregate_fun)
#   return(QF_aggregated)
# }
#
# aggregation_parse <- function(Peptide,
#                               ecols,
#                               fnames,
#                               mode = c("char", "delim"),
#                               chars = 1, sep = "_",
#                               NAnum = NULL,
#                               thr = NULL,
#                               normalize_method = c("quantiles","quantiles.robust","vsn","center.mean",
#                                                    "center.median","div.mean","div.median", "diff.median"),
#                               aggregate_fun = MsCoreUtils::robustSummary){
#   mode <- match.arg(mode)
#   if (any(!apply(Peptide[, ecols], 2, is.numeric))) {
#     stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
#          call. = FALSE)
#   }
#
#   if (tibble::is_tibble(Peptide))
#     Peptide <- as.data.frame(Peptide)
#
#   if (mode == "char") {
#     expdesign <- data.frame(label = colnames(Peptide)[ecols], stringsAsFactors = FALSE) %>%
#       mutate(condition = substr(label, 1, nchar(label) -
#                                   chars), replicate = substr(label, nchar(label) +
#                                                                1 - chars, nchar(label))) %>% unite(ID, condition,
#                                                                                                    replicate, remove = FALSE)
#   }
#   if (mode == "delim") {
#     # colnames(raw) = gsub(get_suffix(colnames(raw)),"", colnames(raw))
#     expdesign <- data.frame(label = colnames(Peptide)[ecols], stringsAsFactors = FALSE) %>%
#       separate(label, c("condition", "replicate"), sep = sep,
#                remove = FALSE, extra = "merge") %>% unite(ID,
#                                                           condition, replicate, remove = FALSE)
#   }
#
#   rownames(expdesign) = expdesign$label
#   Peptide[,ecols] = apply(Peptide[,ecols], 2, function(x){
#     # temp1 <<- x
#     # temp <<- as.integer(x)
#     x[!(!grepl("[A-z]",x) & grepl("\\d",x))] = 0
#     return(as.numeric(x))
#   })
#
#   QF <- readQFeatures(table = Peptide, ecol = ecols,
#                       fnames = fnames,
#                       name = "peptideRaw")
#   rowData(QF[["peptideRaw"]])$nNonZero <- rowSums(assay(QF[["peptideRaw"]]) > 0)
#   if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))
#
#   QF[["peptideRaw"]]@colData = as(expdesign, "DataFrame")
#   colData(QF) = as(expdesign, "DataFrame")
#   QF <- zeroIsNA(QF, "peptideRaw")
#
#
#   filter_thr = function(se,thr){
#     if (any(!c("label", "condition", "replicate") %in%
#             colnames(colData(se)))) {
#       stop("'label', 'condition' and/or 'replicate' columns are not present in '",
#            deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
#            call. = FALSE)
#     }
#     max_repl <- max(colData(se)$replicate)
#     if (thr < 0 | thr > max_repl) {
#       stop("invalid filter threshold applied", "\nRun filter_missval() with a threshold ranging from 0 to ",
#            max_repl)
#     }
#     bin_data <- assay(se)
#     idx <- (is.na(assay(se)))
#     bin_data[!idx] <- 1
#     bin_data[idx] <- 0
#     keep <- bin_data %>% data.frame() %>% rownames_to_column() %>%
#       gather(ID, value, -rowname) %>% left_join(., data.frame(colData(se)),
#                                                 by = "ID") %>% group_by(rowname, condition) %>%
#       summarize(miss_val = n() - sum(value)) %>% filter(miss_val <=
#                                                           thr) %>% spread(condition, miss_val)
#     se_fltrd <- se[keep$rowname, ]
#     return(se_fltrd)
#   }
#
#   if(!is.null(thr))
#     QF[["peptideRaw"]] = filter_thr(QF[["peptideRaw"]], thr = thr)
#   QF <- logTransform(QF, base = 2, i = "peptideRaw", name = "peptideLog")
#   QF <- QFeatures::normalize(QF,
#                              i = "peptideLog",
#                              name = "peptideNorm",
#                              method = normalize_method)
#   # QF[["peptideNorm"]]
#   limma::plotDensities(assay(QF[["peptideLog"]]))
#   limma::plotDensities(assay(QF[["peptideNorm"]]))
#   QF <- filterFeatures(QF, ~ Proteins %in% smallestUniqueGroups(rowData(QF[["peptideLog"]])$Proteins))
#   QF_aggregated <- aggregateFeatures(QF,
#                                      i = "peptideNorm", fcol = "Proteins",
#                                      name = "protein",
#                                      fun = aggregate_fun)
#   return(QF_aggregated)
# }


## extract smallest proteingroups with at least peptide, and exclude the shared peptides with longer protein
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
        sel[i] <- !any(b[[i]] %in% a)#length是1的蛋白名,即included不在list中length > 1的蛋白名中为TRUE，在的话为FALSE
      }
      b <- b[sel] #把length是1的蛋白名,即included不在list中length > 1的蛋白名中为TRUE，提出来，重新赋予b
      j <- j + 1 #j=2再返回while循环
    }
  }

  included <- unlist(included)
  return(included)
}

# make_pe <- function(Peptide, ecols, fnames,
#                     expdesign = NULL,
#                     NAnum = NULL,
#                     thr = NULL){
#   if (any(!apply(Peptide[, ecols], 2, is.numeric))) {
#     stop("specified 'columns' should be numeric", "\nRun make_pe_parse() with the appropriate columns as argument",
#          call. = FALSE)
#   }
#
#   if (any(!c("label", "condition", "replicate") %in%
#           colnames(expdesign))) {
#     stop("'label', 'condition' and/or 'replicate' columns",
#          "are not present in the experimental design",
#          call. = FALSE)
#   }
#
#   if (tibble::is_tibble(Peptide))
#     Peptide <- as.data.frame(Peptide)
#   if (tibble::is_tibble(expdesign))
#     expdesign <- as.data.frame(expdesign)
#
#   rownames(expdesign) = expdesign$label
#
#   raw = Peptide[,ecols]
#   expdesign <- mutate(expdesign, label = make.names(label), condition = make.names(condition)) %>%
#     unite(ID, condition, replicate, remove = FALSE)
#   rownames(expdesign) <- expdesign$ID
#   matched <- match(make.names(delete_prefix(expdesign$label)),
#                    make.names(delete_prefix(colnames(raw))) )
#   if(any(is.na(matched))
#      # && sum(matched,na.rm = T) != length(ecols)
#      ){
#     matched <- match(make.names(expdesign$label),
#                      make.names(colnames(raw)) )
#   }
#   if (any(is.na(matched))
#       # && sum(matched,na.rm = T) != length(ecols)
#       ) {
#     stop("Labels of the experimental design do not match ",
#          "with column names in 'proteins_unique'", "\nRun make_se() with the correct labels in the experimental design",
#          "and/or correct columns specification")
#   }
#   colnames(raw)[matched] <- expdesign$ID
#   raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
#
#   Peptide[,ecols] = apply(raw, 2, function(x){
#     # temp1 <<- x
#     # temp <<- as.integer(x)
#     x[!(!grepl("[A-z]",x) & grepl("\\d",x))] = 0
#     return(as.numeric(x))
#   })
#   colnames(Peptide)[ecols] = colnames(raw)
#   QF <- readQFeatures(table = Peptide, ecol = ecols,
#                       fnames = fnames,
#                       name = "peptideRaw")
#   rowData(QF[["peptideRaw"]])$nNonZero <- rowSums(assay(QF[["peptideRaw"]]) > 0)
#   # if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))
#
#   QF[["peptideRaw"]]@colData = as(expdesign, "DataFrame")
#   colData(QF) = as(expdesign, "DataFrame")
#   QF <- zeroIsNA(QF, "peptideRaw")
#   return(QF)
# }

# make_pe_parse <- function(Peptide,
#                           ecols,
#                           fnames,
#                           mode = c("char", "delim"),
  #                         chars = 1, sep = "_", remove_prefix = T){
  # mode <- match.arg(mode)
  # if (any(!apply(Peptide[, ecols], 2, is.numeric))) {
  #   stop("specified 'columns' should be numeric", "\nRun make_pe_parse() with the appropriate columns as argument",
  #        call. = FALSE)
  # }
  #
  # if (tibble::is_tibble(Peptide))
  #   Peptide <- as.data.frame(Peptide)
  #
  # if(remove_prefix){
  #   colnames(Peptide)[ecols] <- DEP2:::delete_prefix(colnames(Peptide)[ecols]) %>% make.names()
  # }
  # # colnames(Peptide)[ecols] <- DEP2:::delete_prefix(colnames(Peptide)[ecols]) %>% make.names()
  # if (mode == "char") {
  #   expdesign <- data.frame(label = colnames(Peptide)[ecols], stringsAsFactors = FALSE) %>%
  #     mutate(condition = substr(label, 1, nchar(label) -
  #                                 chars), replicate = substr(label, nchar(label) +
  #                                                              1 - chars, nchar(label))) %>% unite(ID, condition,
  #                                                                                                  replicate, remove = FALSE)
  # }
  # if (mode == "delim") {
  #   # colnames(raw) = gsub(get_suffix(colnames(raw)),"", colnames(raw))
  #   expdesign <- data.frame(label = colnames(Peptide)[ecols], stringsAsFactors = FALSE) %>%
  #     separate(label, c("condition", "replicate"), sep = sep,
  #              remove = FALSE, extra = "merge") %>% unite(ID,
  #                                                         condition, replicate, remove = FALSE)
  # }
  #
  # rownames(expdesign) = expdesign$label
  # Peptide[,ecols] = apply(Peptide[,ecols], 2, function(x){
  #   # temp1 <<- x
  #   # temp <<- as.integer(x)
  #   x[!(!grepl("[A-z]",x) & grepl("\\d",x))] = 0
  #   return(as.numeric(x))
  # })
  #
  # QF <- readQFeatures(table = Peptide, ecol = ecols,
  #                     fnames = fnames,
  #                     name = "peptideRaw")
  # rowData(QF[["peptideRaw"]])$nNonZero <- rowSums(assay(QF[["peptideRaw"]]) > 0)
  # # if(!is.null(NAnum))   QF <- filterFeatures(QF, ~ nNonZero >= (ncol(assay(QF[["peptideRaw"]])) - NAnum))
  #
  # QF[["peptideRaw"]]@colData = as(expdesign, "DataFrame")
  # colData(QF) = as(expdesign, "DataFrame")
  # QF <- zeroIsNA(QF, "peptideRaw")
# }

# filter_pe <- function(pe,
#                       thr = NULL,
#                       NAnum = NULL,
#                       filter_column_names = NULL){
#   assertthat::assert_that(class(pe) == "QFeatures", is.null(thr)||is.numeric(thr), is.null(NAnum)||is.numeric(NAnum),
#                           is.null(filter_column_names)||is.character(filter_column_names))
#
#   rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
#   if(!is.null(NAnum))   pe <- filterFeatures(pe, ~ nNonZero >= (ncol(assay(pe[["peptideRaw"]])) - NAnum))
#
#   if(!is.null(thr)){
#     filter_thr = function(se,thr){
#       if (any(!c("label", "condition", "replicate") %in%
#               colnames(colData(se)))) {
#         stop("'label', 'condition' and/or 'replicate' columns are not present in '",
#              deparse(substitute(se)), "'\nRun make_pe() or make_pe_parse() to obtain the required columns",
#              call. = FALSE)
#       }
#       max_repl <- max(colData(se)$replicate)
#       if (thr < 0 | thr > max_repl) {
#         stop("invalid filter threshold applied", "\nRun filter_missval() with a threshold ranging from 0 to ",
#              max_repl)
#       }
#       bin_data <- assay(se)
#       idx <- (is.na(assay(se)))
#       bin_data[!idx] <- 1
#       bin_data[idx] <- 0
#       keep <- bin_data %>% data.frame() %>% rownames_to_column() %>%
#         gather(ID, value, -rowname) %>% left_join(., data.frame(colData(se)),
#                                                   by = "ID") %>% group_by(rowname, condition) %>%
#         summarize(miss_val = n() - sum(value)) %>% filter(miss_val <=
#                                                             thr) %>% spread(condition, miss_val)
#       se_fltrd <- se[keep$rowname, ]
#       return(se_fltrd)
#     }
#     # pe[["peptideRaw"]] = filter_thr(pe[["peptideRaw"]], thr = thr)
#     pe = filter_thr(pe, thr = thr)
#   }
#
#   if(!is.null(filter_column_names)){
#     for(i in filter_column_names){
#       pe <- filterFeatures(pe, as.formula(paste0("~ ", i ,"== ''")))
#     }
#   }
#
#   return(pe)
# }

## Distinguish smallest unique proteingroups and respective unique peptides. Distribute razor peptides to proteingroups
# Peptide_distribution <- function(pe_norm, i = "peptideNorm", fcol = "Proteins"){
#   rd <- rowData(pe_norm[[i]]) %>% as.data.frame()
#   allpgs <- rd[,fcol]
#   SU_PGs <- smallestUniqueGroups(allpgs) ## 得到最终的 smallestUnique PGs
#   SU_PGs_uniquecounts <- allpgs[allpgs %in% SU_PGs] %>% table ## 计算每个最终 PG 的unique数量
#   SU_PGs_uniquecounts %>% class
#   SU_PGs_list1 <- strsplit(SU_PGs,";",fixed = F)
#   names(SU_PGs_list1) <- SU_PGs
#
#   rd$peptide.type <- ifelse(rd[,fcol] %in% SU_PGs, "unique", "razor")
#   nonSU_PGs <- rd %>% filter(peptide.type == "razor") %>% .[,fcol] %>% unique() ## 非 smallestUnique PGs
#   # nonSU_PGs %>% sapply(.,function(x){strsplit(x,";")[[1]]})
#   # nonSU_PGs %>% purrr::map(.,function(x){tem <- strsplit(x,";")[[1]]})
#
#   nonSU_PGs_list <- strsplit(nonSU_PGs,";",fixed = F)
#   names(nonSU_PGs_list) = nonSU_PGs
#
#   # nonSU_PGs_list2 <- nonSU_PGs_list %>% unlist %>% unique()
#
#   (nonSU_PGs_list3 <- nonSU_PGs_list %>% sapply(., function(x){
#     su_pg <- sapply(x, function(y){
#       map(SU_PGs_list1, function(z){
#         if(y %in% z) return(T)
#       }) %>% unlist %>% names
#     }) %>% unlist
#     su_pg_pepcounts <- SU_PGs_uniquecounts[match(su_pg,names(SU_PGs_uniquecounts))]
#     su_pg_pepcounts %>% sort(decreasing = T) %>% names %>% .[1]
#     # %>% names(.) %>% .[1]
#   })) %>% system.time()  ## 这步可能花费大量时间
#   class(nonSU_PGs_list3)
#
#   rd$smallestProteingroups <- ifelse(rd$peptide.type=="unique" , rd[,fcol], nonSU_PGs_list3[match(rd[,fcol],names(nonSU_PGs_list3))])
#   rowData(pe_norm[[i]]) <- rd
#   return(pe_norm)
# }
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




## modify the aggregateQFeatures function in QFeatures, increasing the reserve parameter to reserve certain columns in reduceDataFrame step
.aggregateQFeatures <- function(object, fcol, fun, reserve,...) {
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

  aggregated_assay <- aggregate_by_vector(m, groupBy, fun, ...)
  aggcount_assay <- aggregate_by_vector(m, groupBy, colCounts)
  aggregated_rowdata <- reduceDataFrame(rd, rd[[fcol]],
                                                   simplify = T,
                                                   drop = T,
                                                   count = TRUE,
                                                   reserve = reserve)

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
environment(.aggregateQFeatures) = asNamespace("QFeatures")
assignInNamespace(".aggregateQFeatures", .aggregateQFeatures, ns = "QFeatures")

reduceDataFrame <- function (x, k, count = FALSE, simplify = TRUE, drop = FALSE, reserve)
{
  res <- split(x, k)
  lens <- unname(lengths(res))
  colnames(x)
  if (simplify | drop)
    invars <- QFeatures:::invariant_cols2(res)
  res <- DataFrame(res)

  invars <- c(invars,which(colnames(res) %in% reserve)) %>% unique
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
environment(reduceDataFrame) = asNamespace("QFeatures")
# assignInNamespace("reduceDataFrame", reduceDataFrame, ns = "QFeatures")




####################### for DEP-PTM
# /**
#  * [description]
#  * @param  {[type]} peptidedata         [data frame]
#  * @param  {String} name                [gene name]
#  * @param  {String} ID                  [protein id]
#  * @param  {String} aa                  [Amino acid]
#  * @param  {String} pos                 [Amino acid position in peptide]
#  * @param  {[type]} filter_column_names [the name of filter column]
#  * @param  {[type]} delim               [Sets the delimiter separating the feature names within one protein group]
#  * @return {[type]}                     [description]
#  */
filter_peptide <- function(PTMdata, gene.name = "Gene.names", protein.ID = "Protein",
                           aa = "Amino.acid", pos = "Position",
                           filter_column_names = c("Potential.contaminant", "Reverse"),
                           delim = ";",
                           cutoff_based = "Localization.prob", cutoff = 0.7){
  data = PTMdata
  cols_filt <- match(filter_column_names, colnames(data))

  if (any(is.na(cols_filt))) {
    cols_filt <- cols_filt[!is.na(cols_filt)]
  }
  message("Filtering based on '", paste(filter_column_names,
                                        collapse = "', '"), "' column(s)")
  if (any(!c(gene.name, protein.ID) %in% colnames(data))) {
    stop("'gene.name' and/or 'protein.ID' columns are not present in '",
         deparse(substitute(data)),
         call. = FALSE)
  }

  if (!is.null(cols_filt)) {
    NAs <- is.na(data[, cols_filt])
    data[, cols_filt][NAs] <- ""
    if (length(cols_filt) == 1) {
      data <- dplyr::filter(data, data[, cols_filt] !=
                              "+")
    }
    else if (length(cols_filt) > 1) {
      data <- dplyr::filter(data, !apply(data[, cols_filt] ==
                                           "+", 1, any))
    }
  }


  data$gene_name = ifelse(data[, gene.name] == "" , data[ , protein.ID] , data[ , gene.name]) %>%
    gsub(paste0(delim, ".*"), "", .)
  data$protein_ID = ifelse(data[, protein.ID] == "" , data[ , gene.name] , data[ , protein.ID]) %>%
    gsub(paste0(delim, ".*"), "", .)

  # data %<>% filter(Reverse=="") %>% filter(Potential.contaminant=="")
  data$PTM_ID = paste(data[ , "protein_ID"], "_", data[ , aa] , data[ , pos] ,sep="")
  data$PTM_name = paste(data[, "gene_name"], "_" , data[ , aa] , data[ , pos],sep="")

  if(!is.null(cutoff_based) & cutoff_based != "") {
    cat(paste0("filter based on ", cutoff_based, ": ", cutoff))
    data <- data[data[ , cutoff_based] >= cutoff, ]
  }
  data2 = make_unique(data,names= "PTM_name",ids = "PTM_ID")
  return(data2)
}


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

correct_PTM_by_Protein = function(enriched_peptide , relative_protein, correct_key,
                                  correct_level, unidentified_treatment) {
  assertthat::assert_that(inherits(enriched_peptide, "SummarizedExperiment"),
                          inherits(relative_protein, "SummarizedExperiment"),
                          is.character(correct_key),  length(correct_key) == 1,
                          is.character(unidentified_treatment),  length(unidentified_treatment) == 1)
  if(!correct_key %in% c("name", "ID")){
    stop("correct key should be one of name/ID",call. = F)
  }
  if(!unidentified_treatment %in% c("retain", "remove")){
    stop("unidentified_treatment should be one of retain/remove",call. = F)
  }

  protein_assay = assay(relative_protein)
  protein_df_wide = get_df_wide(relative_protein)
  protein_design <- relative_protein@colData

  key_col = ifelse(correct_key == "name", "gene_name", "protein_ID")

  enriched_peptide2 <<- enriched_peptide
  if(unidentified_treatment == "remove")
    enriched_peptide2 = enriched_peptide2[enriched_peptide2@elementMetadata[, key_col] %in% protein_df_wide[, correct_key],]

  PTM_assay = assay(enriched_peptide2)
  PTM_design <- enriched_peptide2@colData

  protein_assay2 = protein_assay - rowMeans(protein_assay) ## center proteingroup quantity
  protein_assay2 = data.frame(rowname = protein_df_wide[, correct_key] , protein_assay2)
  colnames(protein_assay2)[1] = "key"
  protein_assay2 <- gather(protein_assay2, key = label, value = "pro_expression", -"key")
  protein_assay2[,correct_level] = protein_design[match(protein_assay2$label, protein_design$label),correct_level]
  # protein_assay2$condition = protein_design$condition[match(protein_assay2$label, protein_design$label)]

  ## calculate mean centered quantity of each condition
  thegroups <- c(correct_level, "key")
  protein_assay2_save <<- protein_assay2
  protein_assay3 = protein_assay2 %>%
    group_by_at(vars(one_of(thegroups))) %>%
    # group_by(.dots = thegroups) %>%
    summarise(pro_expression = mean(pro_expression))
  protein_assay3_save <<- protein_assay3

  PTM_assay2 = PTM_assay %>% data.frame() %>% rowid_to_column()
  PTM_assay2[,"key"] = enriched_peptide2@elementMetadata[, key_col]
  PTM_assay2 <- gather(PTM_assay2, key = label, value = "expression", -c("rowid","key"))
  PTM_assay2[,correct_level] = PTM_design[match(PTM_assay2$label, PTM_design$label), correct_level]
  # PTM_assay2$condition = PTM_design$condition[match(PTM_assay2$label, PTM_design$label)]
  PTM_assay2 %>% head()
  PTM_assay2_save <<- PTM_assay2

  ## PTM quantity minus relative proteingroup quantity in different condition
  PTM_assay3 = left_join(PTM_assay2,protein_assay3 , by= c("key",correct_level))
  PTM_assay3$pro_expression[is.na(PTM_assay3$pro_expression)] = 0
  PTM_assay3$expression = PTM_assay3$expression - PTM_assay3$pro_expression
  PTM_assay3_save <<- PTM_assay3
  PTM_assay4 = spread(PTM_assay3[,c("rowid","label","expression")],key = label,value = expression  ) %>%
    arrange(rowid) %>% .[,colnames(enriched_peptide2)] %>% `rownames<-` (rownames(enriched_peptide2))
  PTM_assay4_save <<- PTM_assay4

  assay(enriched_peptide2) = PTM_assay4

  return(enriched_peptide2)

}


# Timecourse_heatmap <- function(ht_mat, exp_design, groupby = "condition",
#                                columns_order = NULL, algo = "cm", k = 6,
#                                color, col_limit, row_font_size,
#                                col_font_size,heatmap_width, heatmap_height,
#                                seed = 10086){
#
#   print("start timeclust")
#   if(is.null(columns_order)) columns_order = exp_design[,groupby] %>% unique()
#
#   # ht_mat = ht_mat_save
#   ht_mat = ht_mat[,which(exp_design[,groupby] %in% columns_order)]
#   ht_mat = ht_mat - rowMeans(ht_mat)
#
#   exp_design2 = exp_design[which(exp_design[,groupby] %in% columns_order),]
#   ht_mat2 <- columns_order %>% map_df(., function(x){
#     rowMeans(ht_mat[,which(exp_design2[,groupby] == x)])
#   }) %>%  t()
#
#   colnames(ht_mat2) = columns_order
#   rownames(ht_mat2) = rownames(ht_mat)
#
#   col_order <- columns_order %>% sapply(., function(x){which(exp_design[,groupby] %in% x)}) %>% unlist
#   col_order_save <<- col_order
#   ht_mat = ht_mat[,col_order]
#
#   set.seed(seed)
#   tca <- timeclust(ht_mat2, algo = algo, k = k, standardize = TRUE)
#   tca_save <<- tca
#   ht_mat2_save <<- ht_mat2
#   ht_mat_save <<- ht_mat
#   cluster1 <- clustCluster(tca)
#   tcadata <- clustData(tca)
#
#   Membershipdata <- clustMembership(tca)
#   mebership <- vector()
#   for (i in 1:length(cluster1)) {
#     mebership[i] <- Membershipdata[names(cluster1)[i],cluster1[i]]
#   }
#   names(mebership) <- names(cluster1)
#
#   res_df <- cbind(name = rownames(ht_mat), Timecourse_cluster = cluster1, mebership = mebership, ht_mat) %>% as.data.frame()
#
#   col_fun_line = circlize::colorRamp2(
#     c(0.01, 0.5,0.9),
#     c("#54c556","#d3ca48","#cd3634"))
#
#   col_fun_heatmap = circlize::colorRamp2(
#     seq(-col_limit, col_limit, (col_limit/5)),
#     rev(RColorBrewer::brewer.pal(11, color)))
#
#   ngroup = length(columns_order)
#
#   panel_fun2 = function(index, nm) {
#     pushViewport(viewport(xscale = c(0,(2*ngroup+1)), yscale = c(-2,2)))
#     grid.rect()
#
#     gby = annotation_axis_grob(at = c(-1.7,0,1.7), labels = c(-1.7,0,1.7), labels_rot = 0,
#                                side = "left", facing = "inside",gp = gpar(fontsize=5))
#     # gby2 = annotation_axis_grob(at = 0,labels = "expression", labels_rot = 90,
#     #                            side = "left", facing = "outside",gp = gpar(fontsize=5))
#     grid.polyline(y=as.vector(t( (tcadata[index,]/4)+0.5 )),
#                   x=c(rep(seq(from=2/(2*ngroup+2),to=(2*ngroup)/(2*ngroup+1),by=2/(2*ngroup+1) ), length(index))),
#                   id=rep(1:length(index), times=rep(ncol(tcadata),times=length(index))),
#                   gp=gpar(col=col_fun_line(mebership[index]), lwd=1,alpha = 0.5))
#
#     if(cluster1[index[1]]==max(cluster1)){
#       gbx = annotation_axis_grob(at = seq(from=2,to=(2*ngroup), by=2), labels = colnames(tcadata), labels_rot = 45,
#                                  side = "bottom", facing = "outside",gp = gpar(fontsize= col_font_size ))
#       grid.draw(gbx)
#     }
#
#     grid.draw(gby)
#     # grid.draw(gby2)
#     popViewport()
#   }
#
#   anno = anno_link(align_to = cluster1, which = "row", panel_fun = panel_fun2, side ="right",
#                    size = unit((heatmap_height * 5/k) - 0.25, "cm"), gap = unit(0, "cm"), width = unit(heatmap_width*1.7, "cm"))
#
#   ht <- ComplexHeatmap::Heatmap(ht_mat,
#                                 heatmap_width = unit(heatmap_width * 5, "cm"),heatmap_height = unit(heatmap_height * 5, "cm"),
#                                 # width = unit(5.5, "cm"),height = unit(20, "cm"),
#                                 col = col_fun_heatmap,
#                                 split = cluster1,
#                                 cluster_rows = F,
#                                 show_row_names = T,
#                                 row_names_side = "left",
#                                 row_dend_width = unit(10, "mm"),
#                                 cluster_row_slices = T,
#                                 cluster_column_slices =F,
#                                 row_gap = unit(1, "mm"),
#                                 #row_km = km,
#                                 cluster_columns =FALSE,
#                                 # column_labels =extractcolname(data),
#                                 heatmap_legend_param = list(
#                                                             title_position = "leftcenter-rot"
#                                                             # height = unit(500,"mm"),
#                                 ),
#                                 row_names_gp = gpar(fontsize = row_font_size),
#                                 column_names_gp = gpar(fontsize = col_font_size),
#                                 right_annotation = rowAnnotation(foo = anno)
#   )
#
#   print("finish heatmap")
#   return(res_save <<- list(res = res_df, ht = ht, mat = ht_mat, timeclust_res = tca))
# }
