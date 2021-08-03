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

#### ===== limma BH FDR ===== #####
test_limma <- function(se, type = c("control", "all", "manual"),
                       control = NULL, test = NULL,
                       design_formula = formula(~ 0 + condition),
                       paired = FALSE) {
  library(limma)
  #require("dplyr", "tidyr", "purrr")

  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  if (paired == FALSE){
    design_formula <- design_formula
  }else{
    design_formula<-formula(~ 0 + condition + replicate)
  }


  # Show error if inputs do not contain required columns
  type <- match.arg(type)

  col_data <- colData(se)
  raw <- assay(se)

  if(any(!c("name", "ID") %in% colnames(rowData(se)))) {
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

  eB_fit <- eBayes(contrast_fit)

  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", adjust.method="BH", coef = comp,
                    number = Inf, confint = TRUE)
    # res <- res[!is.na(res$t),]
    #fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    # res$qval <- res$adj.P.Value
    #res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }

  #limma_res<- topTable(eB_fit, sort.by = 'B', adjust.method="BH", coef = cntrst, number = Inf, confint = T )
  # limma_res$comparison <- rep(cntrst, dim(limma_res)[1])
  #limma_res <- rownames_to_column(limma_res)
  # Retrieve the differential expression test results
  limma_res <- purrr::map_df(cntrst, retrieve_fun)

  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, comparison) %>%
    dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    tidyr::gather(variable, value, -c(rowname,comparison)) %>%
    dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj")) %>%
    tidyr::unite(temp, comparison, variable) %>%
    tidyr::spread(temp, value)
  rowData(se) <- merge(rowData(se), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE)
  return(se)
  #return(table)
}
#manual: logical, if manual heatmap, default: FALSE
#manual_name: act when manual is TRUE, can be one or more contrast format "XXX_vs_XXX", default NULL
#same_trend: same_trend:if the heatmap is the same trend, default TRUE; if the number of maual_name is 2(one from pro1, one from pro2),then the heatmap show: both significant and both trend(the two _vs_), if the number of maual_name is 4,6,8,10,...,then the heatmap show: same trend of significant(any significant), and same trend diff(log2(fc)) from the two part(_vs_), eg: two part is 1,1,1-2,2,2, same signicant can be T,T,T-T,T,T, T,F.F-T,F,F, F,T,F-F,T,F;same trend diff can be (---,---； -++，-++， "-" represent negative value, "+"represent positive and 0 value)
#if_mysplit: if order splices by yourself, default: FALSE
#mysplit: default: NULL, only act when if_mysplit is TRUE and kmeans is TRUE, it is numeric of 1:k, the numeric order is the splices order of heatmap(from top to bottom), eg: c(6,4,1,2,3,5)
#if_rowname_color: if set rowname color which Peptides == 1 to purple, Peptides == 2 to blue, default: FALSE
plot_heatmap <- function (dep, type = c("contrast", "centered"), manual = FALSE, manual_name = NULL, same_trend = FALSE,
                          kmeans = FALSE, k = 6, color = "RdBu", col_limit = 6, indicate = NULL, clustering_distance = c("euclidean",
                                                                                                         "maximum", "manhattan", "canberra",
                                                                                                         "binary", "minkowski", "pearson", "spearman",
                                                                                                         "kendall", "gower"), row_font_size = 6, col_font_size = 10,
                          plot = TRUE, if_mysplit = FALSE, mysplit = NULL, if_rowname_color = FALSE, ...)
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
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(type), is.logical(kmeans), is.numeric(k),
                          length(k) == 1, is.numeric(col_limit), length(col_limit) ==
                            1, is.numeric(row_font_size), length(row_font_size) ==
                            1, is.numeric(col_font_size), length(col_font_size) ==
                            1, is.logical(plot), length(plot) == 1)

  # set rowname of heatmap color, which Peptides == 1 to purple, ==2 to blue
  if(if_rowname_color) {
    rowData(dep)$row_name_color = "black"
    rowData(dep)$row_name_color[which(rowData(dep)$Peptides == 1)] = "purple"
    rowData(dep)$row_name_color[which(rowData(dep)$Peptides == 2)] = "blue"
  }

  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  row_data <- rowData(dep, use.names = FALSE)
  col_data <- colData(dep) %>% as.data.frame()



  if (any(!c("label", "condition", "replicate") %in%
          colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(dep)), "'"), call. = FALSE)
  }
  if (length(grep("_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if (!"significant" %in% colnames(row_data)) {
    stop(paste0("'significant' column is not present in '",
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required column."),
         call. = FALSE)
  }
  if (!is.null(indicate) & type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'",
            call. = FALSE)
  }
  if (!is.null(indicate) & type == "centered") {
    ha1 <- get_annotation(dep, indicate)
  } else {
    ha1 <- NULL
  }
  #filtered <- dep[row_data$significant, ]
  ## Filter for significant proteins only
  library(stringr)
  if(manual) {
    #if it is both trend heatmap, a(data frame), b(check if both significant), you can add a if class(a),class(a) == "DFrame", set constrast muti select,if select two or more, heatmap is desplay both trend heatmap, you can changshi zhe zhong fang fa, tongli ratio, ind ye xugai,zai zai xiafang zhao xugai de difang, eg: all(c(1,0) < 0) | all(c(1,0)>=0)
    #a = row_data[ , paste(manual_name, "_significant", sep = "")]
    #b = apply(a, 1, all)
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
    # if(class(my) == "DFrame"){
    #   ind = sort(unique(as.numeric(unlist(sapply(manual_name, simplify = FALSE, function(x){
    #     which(str_detect(x, gsub("_[[:digit:]]*", "", colname)))
    #   })))))
    # } else {
    #   # when this manual_name can only have one element
    #   # #this is not common use: eg: when colnames(assay(dep)) = "Biotin_1" "Biotin_2" "Biotin_3" "H2A_1"    "H2A_2"    "H2A_3"    "uH2A_1"   "uH2A_2"   "uH2A_3" , set manual_name = "uH2A_vs_Biotin"
    #   ind = which(str_detect(manual_name, gsub("_[[:digit:]]*", "", colname)))
    #   #i = strsplit(manual_name, split = "_")[[1]][c(1,4)]# when DEP shiny app select parse from coloumn,this is c(1,4), else it is c(1,3), because parse from cols, treat_1 hui.bian.cheng treat__1
    #   #ind = unlist(sapply(i, function(x){
    #   #    which(str_detect(colname, x))
    #   #    }))

    # }
    ind = which(gsub("_[[:digit:]]*", "", colname) %in% unlist(strsplit(manual_name, split = "_vs_")))
  } else {
    filtered <- dep[row_data$significant, ]
  }


  # Check for missing values
  if(manual){
    if(any(is.na(assay(filtered)[, ind]))) {
      warning("Missing values in '", deparse(substitute(dep)), "'. ",
              "Using clustering_distance = 'gower'",
              call. = FALSE)
      clustering_distance <- "gower"
      obs_NA <- TRUE
    } else {
      obs_NA <- FALSE
    }

  } else{
    if(any(is.na(assay(filtered)))) {
      warning("Missing values in '", deparse(substitute(dep)), "'. ",
              "Using clustering_distance = 'gower'",
              call. = FALSE)
      clustering_distance <- "gower"
      obs_NA <- TRUE
    } else {
      obs_NA <- FALSE
    }

  }
  # if(any(is.na(assay(filtered)))) {
  #   warning("Missing values in '", deparse(substitute(dep)), "'. ",
  #     "Using clustering_distance = 'gower'",
  #     call. = FALSE)
  #   clustering_distance <- "gower"
  #   obs_NA <- TRUE
  # } else {
  #   obs_NA <- FALSE
  # }

  # Get centered intensity values ('centered')
  if(type == "centered") {
    if(manual){
      rowData(filtered)$mean <- rowMeans(assay(filtered)[ , ind], na.rm = TRUE)
      df <- assay(filtered)[ , ind] - rowData(filtered, use.names = FALSE)$mean
    } else{
      rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
      df <- assay(filtered) - rowData(filtered, use.names = FALSE)$mean
    }
  }
  # Get contrast fold changes ('contrast')
  if(type == "contrast") {
    df <- rowData(filtered, use.names = FALSE) %>%
      data.frame() %>%
      column_to_rownames(var = "name") %>%
      dplyr::select(ends_with("_diff"))
    colnames(df) <-
      gsub("_diff", "", colnames(df)) %>%
      gsub("_vs_", " vs ", .)
    df <- as.matrix(df)
    if(manual){
      i = gsub("_vs_", " vs ", manual_name)
      ii = as.data.frame(df[ , i])
      colnames(ii) = i
      df = as.matrix(ii)
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
    if(type == "centered") {
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
    if(type == "contrast") {
      # Order the k-means clusters according to their average fold change
      order <- data.frame(df) %>%
        cbind(df, cluster = df_kmeans$cluster) %>%
        gather(condition, diff, -cluster) %>%
        group_by(cluster) %>%
        summarize(row = mean(diff)) %>%
        arrange(desc(row)) %>%
        pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }

  # if(ncol(df) == 1) {
  #   col_clust = FALSE
  # } else {
  #   col_clust = TRUE
  # }
  # if(nrow(df) == 1) {
  #   row_clust = FALSE
  # } else {
  #   row_clust = TRUE
  # }
  if(clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }

  # Legend info
  legend <- ifelse(type == "contrast",
                   "log2 Fold change",
                   "log2 Centered intensity")

  if(if_mysplit) {
    cluster_row_slices = FALSE
    cluster_column_slices = FALSE
    split = if(kmeans) {factor(df_kmeans$cluster, levels = mysplit)} else {NULL}
  } else {
    cluster_row_slices = TRUE
    cluster_column_slices = TRUE
    split = if(kmeans) {df_kmeans$cluster} else {NULL}
  }

  #set rowname color
  row_name_color = if(if_rowname_color) {rowData(filtered)$row_name_color} else {"black"}

  # Heatmap
  ht1 = ComplexHeatmap::Heatmap(df,
                col = circlize::colorRamp2(
                  seq(-col_limit, col_limit, (col_limit/5)),
                  rev(RColorBrewer::brewer.pal(11, color))),
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
                row_names_gp = gpar(fontsize = row_font_size, col = row_name_color),
                column_names_gp = gpar(fontsize = col_font_size),
                top_annotation = ha1,
                ...)
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
         "Run plot_pca() with one of the following 'pal' options: ",
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
    stop("'lower' and/or 'upper' arguments are not valid\n         Run plot_pca() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
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
my_plot_volcano <- function(dep, contrast, adjusted = FALSE ,labelWay = "all significant",showNum = 20, chooseTolabel = NULL,fontSize = 9, dotsize = 2, same_width = TRUE, plot = TRUE, fcCutoff, linetype = 4, adjpCutoff = 0.05, if_peptide_color = FALSE, label.rectangle = FALSE, stroke = 0.6, down_color = "#0099CC", stable_color = "#BEBEBE", up_color = "#CD8500", peptide_1_color = "#A020F0", peptide_2_color = "#0000FF") {
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

    df <- data.frame(x = row_data[, diff], y = -log10(row_data[,
        p_values]), significant = row_data[, signif], name = row_data$name, Peptides = row_data$Peptides) %>%
        dplyr::filter(!is.na(significant)) %>% arrange(significant)
    name1 <- gsub("_vs_.*", "", contrast)
    name2 <- gsub(".*_vs_", "", contrast)
    xlimit <- ceiling(max(c(abs(min(df$x)), abs(max(df$x)))))


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
                              palette = c("down" = down_color, 'stable' = stable_color, "up" = up_color, "peptide_1" = peptide_1_color, "peptide_2" = peptide_2_color) ,# #CC9933,#FC4E07,lightcoral,#0099CC,#00AFBB
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
    scale_y_continuous(expand = expand_scale(mult = 0, add = c(0, 0.15))) +
      geom_point(aes(fill = change), data = for_label, size = if(labelWay == "selected proteins"){dotsize + 1} else {dotsize}, shape = 21, color = "black", stroke = stroke, position = "identity")

    if(same_width){
        volcanoPlot <- volcanoPlot + scale_x_continuous(limits = c(-xlimit, xlimit))
    }

    if(adjusted){#for adjust pvalue
          volcanoPlot <- volcanoPlot + geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col=stable_color,lwd=0.6)
      }
      else {#for pvalue
          volcanoPlot <- volcanoPlot + geom_hline(yintercept=min(test$y[!test$change == "stable"]),linetype=linetype, col=stable_color,lwd=0.6)
      #     volcanoPlot <- volcanoPlot + geom_hline(yintercept=-log10(adjpCutoff),linetype=linetype, col="grey",lwd=0.6)
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

#indicate can be one , two, three element, the third means facet_wrap
# if_square: if xlim == ylim, default FALSE
plot_pca <- function (dep, x = 1, y = 2, indicate = c("condition",
    "replicate"), label = FALSE, n = 500, point_size = 4,
    label_size = 3, plot = TRUE, if_square = FALSE)
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
        n, " variable proteins"), x = paste0("PC",
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


#### for gene annotation
library(clusterProfiler)
library(dplyr)
library(parallel)
library(PFAM.db)
library(GO.db)
library(AnnotationDbi)
library(reactome.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
# library(org.Rn.eg.db)

# /**
#  * [description]
#  * @param  {[type]} orgDB [eg: org.Hs.eg.db]
#  * @param  {[character]} gene  the gene name, can be SYMBOL, ENSEMBL, UNIPROT,ALIAS name or mixed of them, note that: when both SYMBOL and ALIAS have a  ENTREZID, will select the SYMBOL mapped ENTREZID]
#  * @return {[type]}       [description]
#  */
my_to_entrezid <- function(orgDB = org.Hs.eg.db, gene) {
  ids1 <- try(mapIds(x = orgDB, keys = gene, keytype = "SYMBOL", column = "ENTREZID"), silent = TRUE)
  ids2 <- try(mapIds(x = orgDB, keys = gene, keytype = "ENSEMBL", column = "ENTREZID"), silent = TRUE)
  ids3 <- try(mapIds(x = orgDB, keys = gene, keytype = "UNIPROT", column = "ENTREZID"), silent = TRUE)
  ids4 <- try(mapIds(x = orgDB, keys = gene, keytype = "ALIAS", column = "ENTREZID"), silent = TRUE)

  ids_lis <- list(ids1 = ids1, ids2 = ids2, ids3 = ids3, ids4 = ids4)

  ids_ind <- c(class(ids1) != "try-error", class(ids2) != "try-error", class(ids3) != "try-error", class(ids4) != "try-error")

  my_ids <- ids_lis[ids_ind]
  my_ids1 <- as.data.frame(do.call(cbind, my_ids), stringsAsFactors = F)

  my_ids1$id <- apply(my_ids1, 1, function(i){
    if(all(is.na(i))) {id = NA
    } else {
      id =  i[which(!is.na(i))[1]]
    }
    return(id)
  })
  return(my_ids1)
}


###############################################
########## geneid : entrezid, class as: character
########## keytype: "ENTREZID"
########## genedb： one of mouse:"org.Mm.eg.db"; human:"org.Hs.eg.db"; rat:org.Rn.eg.db: default "org.Hs.eg.db"
Geneannotate <- function(geneid,keytype="ENTREZID",genedb){
  GOdata <- GO.db
  ##ago <- NA
  #ago<- AnnotationDbi::select(genedb, keys=geneid, columns = c("GO","ENTREZID","ENSEMBL", "SYMBOL", "PFAM", "UNIPROT", "GENENAME", "PATH"), keytype = keytype)
  ago<- AnnotationDbi::select(genedb, keys=geneid, columns = c("GO","ENTREZID", "SYMBOL", "PFAM", "GENENAME", "PATH"), keytype = keytype)
  ago <- ago %>% mutate(numb=c(1:nrow(ago)))
  # try(temgo<-AnnotationDbi::select(GOdata,key=ago$GO,columns = c("DEFINITION","TERM"),keytype = "GOID"))
  try(temgo<-AnnotationDbi::select(GOdata,key=ago$GO,columns = c("TERM"),keytype = "GOID"))
  temgo <- temgo %>% mutate(numb=c(1:nrow(ago)))
  bgo <- merge(ago,temgo,by="numb")

  kegg <- ko2name(ko = paste("ko","",bgo$PATH, sep = ""))
  names(kegg)[2] = "kegg"
  names(kegg)[1] = "PATH"
  kegg[ , 1] = gsub("ko", "", kegg[ , 1])
  kegg = kegg[!duplicated(kegg$kegg), ]

  bgo <- dplyr::left_join(x = bgo, y = kegg, by = "PATH")

  react = AnnotationDbi::select(reactome.db, keys = geneid, columns = c("REACTOMEID", "PATHNAME"), keytypes="ENTREZID")
  names(react)[3] = "react"
  bgo <- dplyr::left_join(x = bgo, y = react, by = "ENTREZID")
  #################### test
  pfamAC = as.data.frame(PFAMID)
  names(pfamAC)[1:2] = c("PFAM", "PFAM_ID")

  bgo = dplyr::left_join(x = bgo, y = pfamAC, by = "PFAM")

  IDS <<- unique(bgo$ENTREZID)
  return(bgo)
}
########### x: Geneannotate函数返回对象bgo$ENTREZID
###########bgo： Geneannotate函数返回对象bgo
mergego <- function(x,bgo){
  library(dplyr)
  subgo <- bgo %>% dplyr::filter(ENTREZID==x)
  subgo <- subgo[!duplicated(subgo$GO),]
  GOCC <- subgo%>%dplyr::filter(ONTOLOGY=="CC") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ; ")
  GOBP <- subgo%>%dplyr::filter(ONTOLOGY=="BP") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ;  ")
  GOMF <- subgo%>%dplyr::filter(ONTOLOGY=="MF") %>%dplyr::select("TERM")%>%unlist()%>% paste0(collapse = " ;  ")

  #DECC <- subgo%>%dplyr::filter(ONTOLOGY=="CC") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #DEBP <- subgo%>%dplyr::filter(ONTOLOGY=="BP") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #DEMF <- subgo%>%dplyr::filter(ONTOLOGY=="MF") %>%dplyr::select("DEFINITION")%>%unlist()%>% paste0(collapse = " ;  ")
  #GOs <- c(GOCC,DECC,GOBP,DEBP,GOMF,DEMF)
  GOs <- c(GOCC,GOBP,GOMF)

  subpfam <- bgo %>% dplyr::filter(ENTREZID==x)
  subpfam <- subpfam[!duplicated(subpfam$PFAM_ID),]
  pfam <- subpfam %>% dplyr::select("PFAM_ID") %>% unlist() %>% paste0(collapse = " ; ")


  sub_gene.descri <- bgo %>% dplyr::filter(ENTREZID==x)
  sub_gene.descri <- sub_gene.descri[!duplicated(sub_gene.descri$GENENAME),]
  gene.descri <- sub_gene.descri %>% dplyr::select("GENENAME") %>% unlist() %>% paste0(collapse = " ; ")

  subkegg <- bgo %>% dplyr::filter(ENTREZID==x)
  subkegg <- subkegg[!duplicated(subkegg$kegg),]
  kegg <- subkegg %>% dplyr::select("kegg") %>% unlist() %>% paste0(collapse = " ; ")

  subreact <- bgo %>% dplyr::filter(ENTREZID==x)
  subreact <- subreact[!duplicated(subreact$react),]
  react <- subreact %>% dplyr::select("react") %>% unlist() %>% paste0(collapse = " ; ")

  res = c(GOs, pfam, gene.descri, kegg, react)
  return(res)
}
############ gene: gene name , class as character
sym_to_entr <- function(gene) {
  ids <- bitr(gene, fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = db)
}

the1stname <- function(gene.name){
  names <- strsplit(as.character(gene.name),";")[[1]][1]
  return(names)
}




########## go analysis function for DEP
library(enrichplot)
library(ggplot2)
library(ggthemes)

library(org.Mm.eg.db)
library(org.Hs.eg.db)
# library(org.Rn.eg.db)
library(DOSE)
# library(openxlsx)
library(GO.db)
library(dplyr)
library(clusterProfiler.dplyr)
#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "human", "mouse" or "rat"
#must be form reat <- goAnalysis(), because following use name reat
#df_with_lg2fc: if df_with_lg2fc
goAnalysis <- function(df, df_with_lg2fc = FALSE, organism="human"){
  if(organism == "human"){
    library(org.Hs.eg.db)
    orgDB <<- org.Hs.eg.db
    #kegg_organism <- "hsa"
  }
  if(organism == "mouse"){
    library(org.Mm.eg.db)
    orgDB <<- org.Mm.eg.db
    #kegg_organism <- "mmu"
  }
  if(organism == "rat"){
    library(org.Rn.eg.db)
    orgDB <<- org.Rn.eg.db
    #kegg_organism <- "rno"
  }


  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }

  # names(ids)[1] = "name"
  #  ids <- inner_join(ids, df, by = "name")

  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)

  ids <- inner_join(ids2, df, by = "name")

  if(df_with_lg2fc){
    ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
    de = ids$fc
    names(de) = unlist(ids$ENTREZID)
    de = sort(de, decreasing = T)
  } else {
    ids <- ids[!is.na(ids$ENTREZID), ]
    de = unlist(ids$ENTREZID)
    names(de) = de
  }

  reat_ALL <- clusterProfiler::enrichGO(gene = names(de), OrgDb = orgDB, ont = "ALL",
                 pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

  reat_BP <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "BP")
  reat_BP@ontology = "BP"

  reat_CC <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "CC")
  reat_CC@ontology = "CC"

  reat_MF <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "MF")
  reat_MF@ontology = "MF"

  reat <- list(ALL = reat_ALL, BP = reat_BP, CC = reat_CC, MF = reat_MF, de = de)
  return(reat)
}

# reat: goAnalysis funtion return value
# ont: one of ALL, BP, CC, MF
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# q.cutoff: qvalue cutoff
# simplify: if remove redundancy of enriched GO terms
giveGO_res_and_table <- function(reat, ont = "BP", pCutoff = 0.05, p.adj.cutoff = 0.05, q.cutoff = 0.2, simplify = FALSE){
  res <- reat[[ont]]
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, qvalue < q.cutoff)

  if(simplify) {
    if(ont == "ALL") {
      sig_res@ontology = "BP"
      sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
      } else {
        sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
      }
  }

  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)

  sig_res@result = sig_res@result[ , c(2:ncol(sig_res@result), 1)]

  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])

  return(table)
}

#legend fold change to log2 fold change
Heatplot <-  function (x, showCategory = 30, foldChange = NULL) {
    n <- update_n(x, showCategory)
    geneSets <- extract_geneSets(x, n)
    foldChange <- fc_readable(x, foldChange)
    d <- list2df(geneSets)
    if (!is.null(foldChange)) {
        d$foldChange <- foldChange[as.character(d[, 2])]
        p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(aes_(fill = ~foldChange),
            color = "white") + scale_fill_continuous(low = "blue",
            high = "red", name = "log2 fold change")
    }
    else {
        p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(color = "white")
    }
    p + xlab(NULL) + ylab(NULL) + theme_minimal() + theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
}
environment(Heatplot) = asNamespace("enrichplot")


Cnetplot <- function (x, showCategory = 5, foldChange = NULL, layout = "kk",
    colorEdge = FALSE, circular = FALSE, node_label = "all",
    ...) {
    node_label <- match.arg(node_label, c("category", "gene",
        "all", "none"))
    if (circular) {
        layout <- "linear"
        geom_edge <- geom_edge_arc
    }
    else {
        geom_edge <- geom_edge_link
    }
    geneSets <- extract_geneSets(x, showCategory)
    g <- list2graph(geneSets)
    foldChange <- fc_readable(x, foldChange)
    size <- sapply(geneSets, length)
    V(g)$size <- min(size)/2
    n <- length(geneSets)
    V(g)$size[1:n] <- size
    if (colorEdge) {
        E(g)$category <- rep(names(geneSets), sapply(geneSets,
            length))
        edge_layer <- geom_edge(aes_(color = ~category), alpha = 0.8)
    }
    else {
        edge_layer <- geom_edge(alpha = 0.8, colour = "darkgrey")
    }
    if (!is.null(foldChange)) {
        fc <- foldChange[V(g)$name[(n + 1):length(V(g))]]
        V(g)$color <- NA
        V(g)$color[(n + 1):length(V(g))] <- fc
        palette <- fc_palette(fc)
        p <- ggraph(g, layout = layout, circular = circular) +
            edge_layer + geom_node_point(aes_(color = ~as.numeric(as.character(color)),
            size = ~size)) + scale_color_gradientn(name = "log2 fold change",
            colors = palette, na.value = "#E5C494")
    }
    else {
        V(g)$color <- "#B3B3B3"
        V(g)$color[1:n] <- "#E5C494"
        p <- ggraph(g, layout = layout, circular = circular) +
            edge_layer + geom_node_point(aes_(color = ~I(color),
            size = ~size))
    }
    p <- p + scale_size(range = c(3, 10), breaks = unique(round(seq(min(size),
        max(size), length.out = 4)))) + theme_void()
    if (node_label == "category") {
        p <- p + geom_node_text(aes_(label = ~name), data = p$data[1:n,
            ], repel = TRUE)
    }
    else if (node_label == "gene") {
        p <- p + geom_node_text(aes_(label = ~name), data = p$data[-c(1:n),
            ], repel = TRUE)
    }
    else if (node_label == "all") {
        p <- p +
            geom_node_text(aes_(label = ~name), data = p$data, repel = TRUE, max.overlaps = Inf)
            # geom_node_text(aes_(label = ~name), data = p$data[-c(1:n),
            # ], repel = TRUE) +
            # geom_node_text(aes_(label = ~name), data = p$data[1:n,
            # ], repel = TRUE)
    }
    return(p)
}

fc_palette <- function (fc)
{
    if (all(fc > 0, na.rm = TRUE)) {
        palette <- color_palette(c("blue", "red"))
    }
    else if (all(fc < 0, na.rm = TRUE)) {
        palette <- color_palette(c("green", "blue"))
    }
    else {
        palette <- color_palette(c("darkgreen", "#0AFF34",
            "#B3B3B3", "#FF6347", "red"))
    }
    return(palette)
}

environment(Cnetplot) = asNamespace("enrichplot")
environment(fc_palette) = asNamespace("enrichplot")

# emapplot
Emapplot <- function (x, showCategory = 30, color = "p.adjust", layout = "kk",
    ...)
{
    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x)
    y <- as.data.frame(x)
    if (is.numeric(n)) {
        y <- y[1:n, ]
    }
    else {
        y <- y[match(n, y$Description), ]
        n <- length(n)
    }
    if (n == 0) {
        stop("no enriched term found...")
    }
    else if (n == 1) {
        g <- graph.empty(0, directed = FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- y$Description
        V(g)$color <- "red"
        return(ggraph(g) + geom_node_point(color = "red",
            size = 5) + geom_node_text(aes_(label = ~name)))
    }
    else {
        id <- y[, 1]
        geneSets <- geneSets[id]
        n <- nrow(y)
        w <- matrix(NA, nrow = n, ncol = n)
        colnames(w) <- rownames(w) <- y$Description
        for (i in 1:n) {
            for (j in i:n) {
                w[i, j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
        wd <- melt(w)
        wd <- wd[wd[, 1] != wd[, 2], ]
        wd <- wd[!is.na(wd[, 3]), ]
        g <- graph.data.frame(wd[, -3], directed = FALSE)
        E(g)$width = sqrt(wd[, 3] * 5)
        g <- delete.edges(g, E(g)[wd[, 3] < 0.2])
        idx <- unlist(sapply(V(g)$name, function(x) which(x ==
            y$Description)))
        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt
        colVar <- y[idx, color]
        V(g)$color <- colVar
    }
    p <- ggraph(g, layout = layout)
    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha = 0.8, aes_(width = ~I(width)),
            colour = "darkgrey")
    }
    p + geom_node_point(aes_(color = ~color, size = ~size)) +
        geom_node_text(aes_(label = ~name), repel = TRUE) + theme_void() +
        scale_color_continuous(low = "red", high = "blue",
            name = color, guide = guide_colorbar(reverse = TRUE)) +
        scale_size(range = c(3, 8))
}
environment(Emapplot) = asNamespace("enrichplot")

#plot_type: barplot, dotplot,emaplot, cnetplot, plotGOgraph, goplot(based on your significant limit)

# res: giveGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
my_barplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(!Split){
      sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
      barplot(sig_res_for_ALL, showCategory = ShowCategory, color = color) + facet_grid(ONTOLOGY~., scales = "free", space = "free")
      } else {
        barplot(res$sig_res, showCategory = ShowCategory, color = color, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free", space = "free")
      }
    } else {
      barplot(res$sig_res, showCategory = ShowCategory, color = color)
    }
}


# res: giveGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
my_dotplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(!Split){
      sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
      dotplot(sig_res_for_ALL, showCategory = ShowCategory, color = color) + facet_grid(ONTOLOGY~., scales = "free", space = "free")
      } else {
        dotplot(res$sig_res, showCategory = ShowCategory, color = color, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free", space = "free")
      }
    } else {
      dotplot(res$sig_res, showCategory = ShowCategory, color = color)
    }
}

# res: giveGO_res_and_table funtion return value
# color: p.adjust" or "pvalue"
# ShowCategory: number of go terms of plot
my_dotplot_opt <- function(res, color = "p.adjust", size = "Count", title = "", decreasing = TRUE, ShowCategory = 20){
  df = as.data.frame(res$sig_res)

  if(nrow(df) > ShowCategory) {
    df = df[1:ShowCategory, ]
  }
  df$x = ""
  orderBy = color
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels = unique(df$Description[idx]))

p_up = as.numeric(format(max(df[ , color]), digits = 2))
p_low = as.numeric(format(min(df[ , color]), digits = 2))

if(color == "p.adjust"){
  col.name = "FDR"
}
if(color == "pvalue"){
  col.name = "Pvalue"
}


a =  ggplot(df, aes_string(x = "x", y = "Description", size = size, color = color)) +
  geom_point() +
  geom_point(shape = 21, color = "black", stroke = 0.6, position = "identity")+
  scale_colour_gradient2(low = "#df293f", mid = ifelse(!max(df[ , color]) == min(df[ , color]), "#f9f9f9", "#f19c9b"), high = "#959cc8", midpoint = min(df[ , color]) + (max(df[ , color])-min(df[ , color]))/2, name = col.name, breaks = unique(c(min(df[ , color]), max(df[ , color]))), limits = c(min(df[ , color]), max(df[ , color])), labels = unique(c(format(p_low, scientific = TRUE), format(p_up, scientific = TRUE)))) +
  ylab(NULL) +
  ggtitle(title) +
  scale_size(name = "Number of proteins found", range = c(3, 8), breaks = c(min(df[ , "Count"]), max(df[ , "Count"])), limits = c(min(df[ , "Count"]), max(df[ , "Count"]))) +
  theme(panel.background=element_blank(),
        axis.ticks=element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key=element_blank(),
        legend.spacing = unit(3, "lines"),
        legend.justification="bottom", legend.position=c(0.15,0),#0.15
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        #plot.margin = margin(t = 0, r = 0, b = 0.1, l = 15, unit = "cm"),
        title = element_text(size = 10)#all title elements: plot, axes, legends
  )+
  labs(x = NULL) +
  guides(
    color = guide_colorbar(order = 1, frame.colour = "black", frame.linetype = 1, frame.linewidth = 1, draw.ulim = TRUE, draw.llim = TRUE, raster = TRUE, ticks = FALSE, label = TRUE, title.position = "left", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(angle = 90,  hjust = 0.5, size = 10), reverse = TRUE, barwidth = unit(0.5, units = "cm"), barheight = unit(2.5, units = "cm")),
    size = guide_legend(order = 0, title.position = "right", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(size = 10), override.aes = list(alpha = 1, bg = "grey", color = "black"), reverse = TRUE, keyheight = 4)
  ) +
  scale_x_discrete(expand = expand_scale(mult = 0, add = c(0.008, 0.25)))+#0.005, 0.25
  scale_y_discrete()
return(a)
}

# res: giveGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
#  ?ggraph
my_emaplot <- function(res, ShowCategory = 30, color = "p.adjust", layout = "kk", ont = "BP") {
  if(ont == "ALL") {
    sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
  try(Emapplot(sig_res_for_ALL, showCategory = ShowCategory, color = color, layout = layout), silent = T)
    } else {
  try(Emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = layout), silent = T)
    }
}


# res: giveGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# circular: if circular layout
# colorEdge: if colorEdge
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
my_cnetplot <- function(res, ShowCategory = 5, circular = TRUE, colorEdge = TRUE, df_with_lg2fc = FALSE, ont = "BP") {
  if(ont == "ALL") {
    sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
      if(df_with_lg2fc) {
    try(Cnetplot(x = sig_res_for_ALL, showCategory = ShowCategory, foldChange = res$de, circular = circular, colorEdge = colorEdge), silent = T)
  } else {
    try(Cnetplot(x = sig_res_for_ALL, showCategory = ShowCategory, foldChange=NULL, circular = circular, colorEdge = colorEdge), silent = T)
  }

    } else {
      if(df_with_lg2fc) {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange = res$de, circular = circular, colorEdge = colorEdge), silent = T)
  } else {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange=NULL, circular = circular, colorEdge = colorEdge), silent = T)
  }
    }
}

# res: giveGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# ont: the same as function giveGO_res_and_table, "ALL","BP","CC", or "MF"
# Layout: layout
# circular: if circular
my_goplot <- function(res, ShowCategory = 10, color = "p.adjust", ont = "BP", Layout = "kk",circular = TRUE) {
  library(ggplotify)
    if(! ont == "ALL"){
    a = try(goplot(res$sig_res, showCategory = ShowCategory,color = color, layout = Layout,geom = "text", circular = circular), silent = T)
    if(!class(a) == "try-error") {
      print(a)
      } else {
        print("can not plot")
      }
  } else {
      print("Ontology ALL: can not plot for goplot")
    }
}



# res: giveGO_res_and_table funtion return value
# firstSigNodes: number of significant nodes (retangle nodes in the graph)
my_plotGOgraph <- function(res, firstSigNodes = 10, ont = "BP") {
  library(ggplotify)
    if(! ont == "ALL"){
    a = function() {try(plotGOgraph(res$sig_res, firstSigNodes = firstSigNodes), silent = T)}
    a = try(as.ggplot(a), silent = T)
    if(!class(a) == "try-error") {
      print(a)
      } else {
        print("can not plot")
      }
  } else {
      print("Ontology ALL: can not plot for plotGOgraph")
    }
}


# res: giveGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
my_heatplot <- function(res, ShowCategory = 30, df_with_lg2fc = FALSE, ont = "BP") {
  if(ont == "ALL") {
    sig_res_for_ALL <- res$sig_res %>% clusterProfiler.dplyr::arrange(p.adjust)
      if(df_with_lg2fc) {
    try(Heatplot(sig_res_for_ALL, showCategory = ShowCategory, foldChange = res$de), silent = T)
    } else {
    try(Heatplot(sig_res_for_ALL, showCategory = ShowCategory, foldChange = NULL), silent = T)
    }
    } else {
        if(df_with_lg2fc) {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)
    } else {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = NULL), silent = T)
    }
    }
}


########## kegg and reactome analysis function for DEP
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
# library(org.Rn.eg.db)
library(DOSE)
# library(openxlsx)
library(GO.db)
library(dplyr)
library(ReactomePA)

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "hsa", "mmu", "rno"
#must be form reat <- keggAnalysis(), because following use name reat
keggAnalysis <- function(df, organism="hsa", df_with_lg2fc = FALSE){
  if(organism == "hsa"){
    library(org.Hs.eg.db)
    orgDB <<- org.Hs.eg.db
    #kegg_organism <- "hsa"
  }
  if(organism == "mmu"){
    library(org.Mm.eg.db)
    orgDB <<- org.Mm.eg.db
    #kegg_organism <- "mmu"
  }
  if(organism == "rno"){
    library(org.Rn.eg.db)
    orgDB <<- org.Rn.eg.db
    #kegg_organism <- "rno"
  }


  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }

  # names(ids)[1] = "name"
  #  ids <- inner_join(ids, df, by = "name")

  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)

  ids <- inner_join(ids2, df, by = "name")

  if(df_with_lg2fc){
    ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
    de = ids$fc
    names(de) = unlist(ids$ENTREZID)
    de = sort(de, decreasing = T)
  } else {
    ids <- ids[!is.na(ids$ENTREZID), ]
    de = unlist(ids$ENTREZID)
    names(de) = de
  }

  reat <- enrichKEGG(gene = names(de), organism = organism,
                 pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = de)
  return(reat)
}

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: one of "human", "rat", "mouse"
#must be form reat <- keggAnalysis(), because following use name reat
reactAnalysis <- function(df, organism="human", df_with_lg2fc = FALSE){
  if(organism == "human"){
    library(org.Hs.eg.db)
    orgDB <<- org.Hs.eg.db
    #kegg_organism <- "hsa"
  }
  if(organism == "mouse"){
    library(org.Mm.eg.db)
    orgDB <<- org.Mm.eg.db
    #kegg_organism <- "mmu"
  }
  if(organism == "rat"){
    library(org.Rn.eg.db)
    orgDB <<- org.Rn.eg.db
    #kegg_organism <- "rno"
  }


  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }

  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")

  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)

  ids <- inner_join(ids2, df, by = "name")

  if(df_with_lg2fc){
    ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
    de = ids$fc
    names(de) = unlist(ids$ENTREZID)
    de = sort(de, decreasing = T)
  } else {
    ids <- ids[!is.na(ids$ENTREZID), ]
    de = unlist(ids$ENTREZID)
    names(de) = de
  }

  reat <- ReactomePA::enrichPathway(gene = names(de), organism = organism,
                 pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable=TRUE)

  reat <- list(res = reat, de = de)
  return(reat)
}

# reat: keggAnalysis or reactAnalysis funtion return value
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# q.cutoff: qvalue cutoff
givekegg_reat_res_and_table <- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.05, q.cutoff = 0.2){
  res <- reat$res
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, qvalue < q.cutoff)

  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)

  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])

  return(table)
}

#plot_type: barplot, dotplot,emaplot, cnetplot(based on your significant limit)
# res_kegg or res_react: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# barplot(res_kegg$sig_res, showCategory = ShowCategory, color = color)
# barplot(res_react$sig_res, showCategory = ShowCategory, color = color)

# res_kegg or res_react: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# dotplot(res_kegg$sig_res, showCategory = ShowCategory, color = color)
# dotplot(res_react$sig_res, showCategory = ShowCategory, color = color)

# res: givekegg_reat_res_and_table or giveGO_res_and_table funtion return value
# color: p.adjust" or "pvalue"
# ShowCategory: number of go terms of plot
my_dotplot_opt <- function(res, color = "p.adjust", size = "Count", title = "", decreasing = TRUE, ShowCategory = 20){
  df = as.data.frame(res$sig_res)

  if(nrow(df) > ShowCategory) {
    df = df[1:ShowCategory, ]
  }
  df$x = ""
  orderBy = color
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels = unique(df$Description[idx]))

p_up = as.numeric(format(max(df[ , color]), digits = 2))
p_low = as.numeric(format(min(df[ , color]), digits = 2))

if(color == "p.adjust"){
  col.name = "FDR"
}
if(color == "pvalue"){
  col.name = "Pvalue"
}


a =  ggplot(df, aes_string(x = "x", y = "Description", size = size, color = color)) +
  geom_point() +
  geom_point(shape = 21, color = "black", stroke = 0.6, position = "identity")+
  scale_colour_gradient2(low = "#df293f", mid = ifelse(!max(df[ , color]) == min(df[ , color]), "#f9f9f9", "#f19c9b"), high = "#959cc8", midpoint = min(df[ , color]) + (max(df[ , color])-min(df[ , color]))/2, name = col.name, breaks = unique(c(min(df[ , color]), max(df[ , color]))), limits = c(min(df[ , color]), max(df[ , color])), labels = unique(c(format(p_low, scientific = TRUE), format(p_up, scientific = TRUE)))) +
  ylab(NULL) +
  ggtitle(title) +
  scale_size(name = "Number of proteins found", range = c(3, 8), breaks = c(min(df[ , "Count"]), max(df[ , "Count"])), limits = c(min(df[ , "Count"]), max(df[ , "Count"]))) +
  theme(panel.background=element_blank(),
        axis.ticks=element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key=element_blank(),
        legend.spacing = unit(3, "lines"),
        legend.justification="bottom", legend.position=c(0.15,0),#0.15
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        #plot.margin = margin(t = 0, r = 0, b = 0.1, l = 15, unit = "cm"),
        title = element_text(size = 10)#all title elements: plot, axes, legends
  )+
  labs(x = NULL) +
  guides(
    color = guide_colorbar(order = 1, frame.colour = "black", frame.linetype = 1, frame.linewidth = 1, draw.ulim = TRUE, draw.llim = TRUE, raster = TRUE, ticks = FALSE, label = TRUE, title.position = "left", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(angle = 90,  hjust = 0.5, size = 10), reverse = TRUE, barwidth = unit(0.5, units = "cm"), barheight = unit(2.5, units = "cm")),
    size = guide_legend(order = 0, title.position = "right", title.theme = element_text(angle = 90,  hjust = 0.5, size = 10), label.theme = element_text(size = 10), override.aes = list(alpha = 1, bg = "grey", color = "black"), reverse = TRUE, keyheight = 4)
  ) +
  scale_x_discrete(expand = expand_scale(mult = 0, add = c(0.008, 0.25)))+#0.005, 0.25
  scale_y_discrete()
return(a)
}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
emaplot_for_react_kegg <- function(res, ShowCategory = 30, color = "p.adjust", layout = "kk") {
  try(Emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = layout), silent = T)
}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# circular: if circular layout
# colorEdge: if colorEdge
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
cnetplot_for_react_kegg <- function(res, ShowCategory = 5, circular = TRUE, colorEdge = TRUE, df_with_lg2fc = FALSE) {
    if(df_with_lg2fc) {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange = res$de, circular = circular, colorEdge = colorEdge), silent = T)
  } else {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange=NULL, circular = circular, colorEdge = colorEdge), silent = T)
  }

}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
heatplot_for_react_kegg <- function(res, ShowCategory = 30, df_with_lg2fc = FALSE) {
    if(df_with_lg2fc) {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)
    } else {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = NULL), silent = T)
    }
}

########## gsego analysis function for DEP
library(enrichplot)
library(ggplot2)
library(ggthemes)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
# library(org.Rn.eg.db)
library(DOSE)
# library(openxlsx)
library(GO.db)
library(dplyr)
library(clusterProfiler.dplyr)
#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "human", "mouse" or "rat"
#must be form reat <- goAnalysis(), because following use name reat
gsegoAnalysis <- function(df, organism="human"){
  if(organism == "human"){
    library(org.Hs.eg.db)
    orgDB <<- org.Hs.eg.db
    #kegg_organism <- "hsa"
  }
  if(organism == "mouse"){
    library(org.Mm.eg.db)
    orgDB <<- org.Mm.eg.db
    #kegg_organism <- "mmu"
  }
  if(organism == "rat"){
    library(org.Rn.eg.db)
    orgDB <<- org.Rn.eg.db
    #kegg_organism <- "rno"
  }


  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }

  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")

  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)

  ids <- inner_join(ids2, df, by = "name")

  ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  de = ids$fc
  names(de) = unlist(ids$ENTREZID)
  de = sort(de, decreasing = T)

  set.seed(1234)
  reat_ALL <- try(gseGO(gene = de, OrgDb = orgDB, ont = "ALL",
                 pAdjustMethod = "BH", nPerm = 1000, pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE)
  reat_ALL <- setReadable(reat_ALL, OrgDb = orgDB, keyType="ENTREZID")

  reat_BP <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "BP")
  reat_BP@setType = "BP"

  reat_CC <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "CC")
  reat_CC@setType = "CC"

  reat_MF <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "MF")
  reat_MF@setType = "MF"

  reat <- list(ALL = reat_ALL, BP = reat_BP, CC = reat_CC, MF = reat_MF, de = de)
  return(reat)
}

# reat: gsegoAnalysis funtion return value
# ont: one of ALL, BP, CC, MF
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# NES.cutoff: NES cutoff,it means abs(NES) > NES.cutoff
# simplify: if remove redundancy of enriched GO terms
# Phenotype: the Phenotype you want to show, one or both of c("activated", "suppressed"), default: c("activated", "suppressed")
givegseGO_res_and_table <- function(reat, ont = "BP", pCutoff = 0.05, p.adj.cutoff = 0.25, NES.cutoff = 1, simplify = FALSE, Phenotype = c("activated", "suppressed")){
  res <- reat[[ont]]
  res@result$phenotype = ifelse(res@result$NES > 0, "activated", "suppressed")
  res@result$group = paste(res@result$ONTOLOGY, res@result$phenotype, sep = "_")

  if(length(Phenotype) == 1) {
    res <- clusterProfiler.dplyr::filter(res, phenotype == Phenotype)#
  }

  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, abs(NES) > NES.cutoff)

  if(simplify) {
    if(ont == "ALL") {
      sig_res@setType = "BP"
      sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
    } else {
      sig_res <- clusterProfiler::simplify(sig_res, cutoff=0.7, by="p.adjust", select_fun=min)
    }
  }

  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)

  sig_res@result = sig_res@result[ , c(2:ncol(sig_res@result), 1)]

  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])

  return(table)
}

#plot_type: barplot, dotplot,emaplot, cnetplot, heatplot, gseaplot2

# res: givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
gse_dotplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(Split) {
      enrichplot:::dotplot(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "group") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free")
    } else {
      enrichplot:::dotplot(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "phenotype") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free")
    }
  } else {
    enrichplot:::dotplot(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split="group")
  }
}

# res: givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
gse_barplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(Split) {
      enrichplot:::barplot.enrichResult(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "group") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free") + labs(y = "NES")
    } else {
      enrichplot:::barplot.enrichResult(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split = "phenotype") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free") + labs(y = "NES")
    }
  } else {
    enrichplot:::barplot.enrichResult(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split="group") + labs(y = "NES")
  }
}

# res: giveGO_res_and_table or givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# try(emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = "kk"), silent = T)


# res: giveGO_res_and_table or givegseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# foldChange: res$de
# circular: TRUE or FALSE
# try(cnetplot(res$sig_res, showCategory = ShowCategory, color = color, foldChange = res$de, circular = circular, colorEdge = TRUE), silent = T)

# ShowCategory: number of go terms of plot
# foldChange: res$de
# try(heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)

#enrichplot:::tableGrob2 修改，将ptable的color去掉,使description的颜色由绿色变为黑色，仅适用于ES有一条线
mytableGrob2 <- function (d, p = NULL)
{
  d <- d[order(rownames(d)), ]
  tp <- gridExtra::tableGrob(d)#, theme = gridExtra::ttheme_default(base_size = 9),theme = gridExtra::ttheme_minimal()
  if (is.null(p)) {
    return(tp)
  }
  pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  j <- which(tp$layout$name == "rowhead-fg")
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i + 1]][["gp"]] = grid::gpar(col = "black", cex = 0.8)
  }
  return(tp)
}


my_gseaplot2 <- function (x, geneSetID, title = "", color = "green", base_size = 11,
          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE,
          ES_geom = "line")
{
  library(ggplot2)
  library(dplyr)
  # library(openxlsx)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- enrichplot:::gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description),
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description),
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8),
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(),
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(),
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2,
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description ==
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin,
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) +
    ylab(NULL) + theme_classic(base_size) + theme(legend.position = "none",
                                                  plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
                                                  axis.ticks = element_blank(), axis.text = element_blank(),
                                                  axis.line.x = element_blank()) + scale_x_continuous(expand = c(0,
                                                                                                                 0)) + scale_y_continuous(expand = c(0, 0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0)
      inv <- inv + 1
    col = c(rev(RColorBrewer::brewer.pal(5, "Blues")), RColorBrewer::brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[unique(inv)])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin,
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax,
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d,
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x,
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked list metric") + xlab("Rank in Ordered Dataset") +
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2,
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust", "NES")]
    pd <- pd[order(pd[, 1], decreasing = FALSE), ]
    rownames(pd) <-  ""

    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- mytableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp,
                                                                         xmin = quantile(p.res$data$x, 0.7), xmax = quantile(p.res$data$x,
                                                                                                                             0.95), ymin = quantile(p.res$data$runningScore,
                                                                                                                                                    0.75), ymax = quantile(p.res$data$runningScore,
                                                                                                                                                                           0.9))
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(),
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2,
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  PloT <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
  PloT + theme(plot.margin = margin(t = 1,
                                     r = 1, b = 0.2, l = 0.2, unit = "cm"))
}

########## gseKEGG gsereactome analysis function for DEP
library(enrichplot)
library(ggplot2)
library(ggthemes)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
# library(org.Rn.eg.db)
library(DOSE)
# library(openxlsx)
library(GO.db)
library(dplyr)
library(clusterProfiler.dplyr)
library(ReactomePA)


#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "hsa", "mmu", "rno"
#must be form reat <- gsekeggAnalysis(), because following use name reat
gsekeggAnalysis <- function(df, organism="hsa"){
  if(organism == "hsa"){
    library(org.Hs.eg.db)
    orgDB <<- org.Hs.eg.db
    #kegg_organism <- "hsa"
  }
  if(organism == "mmu"){
    library(org.Mm.eg.db)
    orgDB <<- org.Mm.eg.db
    #kegg_organism <- "mmu"
  }
  if(organism == "rno"){
    library(org.Rn.eg.db)
    orgDB <<- org.Rn.eg.db
    #kegg_organism <- "rno"
  }


  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }

  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")

  #   de = ids$fc
  #   names(de) = ids$ENTREZID
  #   de = sort(de, decreasing = T)

  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)

  ids <- inner_join(ids2, df, by = "name")

  ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  de = ids$fc
  names(de) = unlist(ids$ENTREZID)
  de = sort(de, decreasing = T)

    set.seed(1234)
  reat <- try(gseKEGG(gene = de, organism = organism,
                 pAdjustMethod = "BH", nPerm = 1000, pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = de)
  return(reat)
}

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: one of "human", "rat", "mouse"
#must be form reat <- gsereactAnalysis(), because following use name reat
gsereactAnalysis <- function(df, organism="human"){
  if(organism == "human"){
    library(org.Hs.eg.db)
    orgDB <<- org.Hs.eg.db
    #kegg_organism <- "hsa"
  }
  if(organism == "mouse"){
    library(org.Mm.eg.db)
    orgDB <<- org.Mm.eg.db
    #kegg_organism <- "mmu"
  }
  if(organism == "rat"){
    library(org.Rn.eg.db)
    orgDB <<- org.Rn.eg.db
    #kegg_organism <- "rno"
  }


  # ids <- bitr(df$name, fromType = "SYMBOL",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB)
  # ids1 <- try(bitr(setdiff(df$name, ids$SYMBOL), fromType = "UNIPROT",
  #             toType = "ENTREZID",
  #             OrgDb = orgDB))
  # if(class(ids1) != "try-error") {
  #   names(ids1)[1] = "SYMBOL"
  #   ids = rbind(ids, ids1)
  # }

  # names(ids)[1] = "name"
  # ids <- inner_join(ids, df, by = "name")

  #   de = ids$fc
  #   names(de) = ids$ENTREZID
  #   de = sort(de, decreasing = T)

  ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)

  ids <- inner_join(ids2, df, by = "name")

  ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  de = ids$fc
  names(de) = unlist(ids$ENTREZID)
  de = sort(de, decreasing = T)

  set.seed(1234)
  reat <- try(ReactomePA::gsePathway(gene = de, organism = organism,
                 pAdjustMethod = "BH", nPerm = 1000, pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = de)
  return(reat)
}


# reat: gsekeggAnalysis or gsereactAnalysis funtion return value
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# NES.cutoff: NES cutoff,it means abs(NES) > NES.cutoff
# Phenotype: the Phenotype you want to show, one or both of c("activated", "suppressed"), default: c("activated", "suppressed")
give_gsekegg_gsereat_res_and_table <- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.25, NES.cutoff = 1, Phenotype = c("activated", "suppressed")){
  res <- reat$res
  res@result$phenotype = ifelse(res@result$NES > 0, "activated", "suppressed")

  if(length(Phenotype) == 1) {
    res <- clusterProfiler.dplyr::filter(res, phenotype == Phenotype)#
  }

  sig_res <- clusterProfiler.dplyr::filter(res, pvalue < pCutoff, p.adjust < p.adj.cutoff, abs(NES) > NES.cutoff)

  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)

  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])

  return(table)
}

#plot_type: barplot, dotplot,emaplot, cnetplot, heatplot, gseaplot2

# res_gsekegg or res_gsereact: give_gsekegg_gsereat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
#enrichplot:::barplot.enrichResult(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split="phenotype") + labs(y = "NES")

# res_gsekegg or res_gsereact: give_gsekegg_gsereat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# enrichplot:::dotplot(res$sig_res, showCategory = ShowCategory, color = color, x = "NES", split="phenotype")

# res: give_gsekegg_gsereat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# try(emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = "kk"), silent = T)


# res: give_gsekegg_gsereat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# foldChange: res$de
# circular: TRUE or FALSE
# try(cnetplot(res$sig_res, showCategory = ShowCategory, color = color, foldChange = res$de, circular = circular, colorEdge = TRUE), silent = T)

# ShowCategory: number of go terms of plot
# foldChange: res$de
# try(heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)

## for RNAseq differential analysis
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(BiocParallel)
library(apeglm)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(IHW)
library(ashr)
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
                "Run plot_pca() with 'x' and 'y' <= ", ncol(data),
                "."), call. = FALSE)
  }
  if (n > nrow(data)) {
    stop(paste0("'n' argument is not valid.\n", "Run plot_pca() with 'n' <= ",
                nrow(data), "."), call. = FALSE)
  }
  # for shiny indicate can select the third indicate
  mycolData <- coldata
  mycolData$Condition = mycolData$condition
  mycolData$Replicate = mycolData$replicate

  columns <- colnames(mycolData)
  if (!is.null(indicate)) {
    if (length(indicate) > 3) {
      stop("Too many features in 'indicate'\n        Run plot_pca() with a maximum of 3 indicate features")
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
  labs(title = paste0("PCA plot - top ", n, " variable proteins"), x = paste0("PC", x, ": ", percent[x], "%"), y = paste0("PC", y, ": ", percent[y], "%")) +
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
         "Run plot_pca() with one of the following 'pal' options: ",
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
    stop("'lower' and/or 'upper' arguments are not valid\n         Run plot_pca() with 'lower' and 'upper' between -1 and 1",
         call. = FALSE)
  }
  pals <- RColorBrewer::brewer.pal.info %>% rownames_to_column() %>%
    filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel",
         " (qualitative panels also not allowed)\n",
         "Run plot_pca() with one of the following 'pal' options: ",
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
plot_heatmap_rnaseq <- function (data, row_data, type = c("raw", "centered"), manual = FALSE, manual_name = NULL, kmeans = FALSE, k = 6, color = "RdBu", col_limit = 6, clustering_distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall", "gower"), row_font_size = 6, col_font_size = 10,plot = TRUE, if_mysplit = FALSE, mysplit = NULL, ...)
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

  # Heatmap
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
  return(res_df)
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
  assertthat::assert_that(is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),length(lfc) == 1)
  row_data <- as.data.frame(diff)

  if (length(grep("_padj|_log2FoldChange", colnames(row_data))) <
      1) {
    stop("'[contrast]_log2FoldChange' and/or '[contrast]_padj' columns are not present in '",
         deparse(substitute(diff)), "'\nRun get_res() to obtain the required columns",
         call. = FALSE)
  }
  cols_p <- grep("_padj", colnames(row_data))
  cols_diff <- grep("_log2FoldChange", colnames(row_data))
  if (length(cols_p) == 1) {
    diff$significant <- row_data[, cols_p] <= alpha &
      abs(row_data[, cols_diff]) >= lfc
    diff$significant[which(is.na(diff$significant))] = FALSE
    diff$contrast_significant <- diff$significant
    colnames(diff)[ncol(diff)] <- gsub("padj","significant", colnames(row_data)[cols_p])
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

    diff <- merge(diff, sign_df, by = "row.names")
    rownames(diff) = diff$Row.names
    diff = diff[,-1]
  }
  # diff$significant[which(is.na(diff$significant))] = FALSE
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
plot_pro_rna_heatmap <- function(rna, pro, gene, Type = "all", row_font_size = 6, col_font_size = 8, column_title_size = 10, Numcol_pro = 2, Numcol_rna = 2, cluster_columns = F, color = "RdBu") {
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
plot_pro_rna_heatmap_for_whole <- function(rna, pro, manual = FALSE, contrast = NULL, Type = "at least one", row_font_size = 6, col_font_size = 8, column_title_size = 10, Numcol_pro = 2, Numcol_rna = 2, cluster_columns = F, color = "RdBu") {
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


## speed up the speed of GO analysis
kegg_list <- function (db) {
  url <- paste0("http://rest.kegg.jp/list/", db, collapse = "")
  clusterProfiler:::kegg_rest(url)
}

ko2name <- function (ko) {
  p <- clusterProfiler:::kegg_list("pathway")
  ko2 <- gsub("^ko", "path:map", ko)
  ko.df <- data.frame(ko = ko, from = ko2)
  res <- merge(ko.df, p, by = "from", all.x = TRUE)
  res <- res[, c("ko", "to")]
  colnames(res) <- c("ko", "name")
  return(res)
}


get_GOdata_Env <- function () {
  if (!exists(".GO_data_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_data_Env", new.env(), envir=envir)
  }
  get(".GO_data_Env", envir = .GlobalEnv)
}

enrichGO <- function(gene,
                     OrgDb,
                     keyType = "ENTREZID",
                     ont="MF",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     universe,
                     qvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 500,
                     readable=FALSE, pool=FALSE) {
  ont %<>% toupper
  ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

  OrgDb_name = OrgDb$packageName
  OrgDb_version = package.version(OrgDb$packageName)
  GO_DATAfolder = paste(system.file(package = "DEP2"),"/GOdata",sep="")
  if (!file.exists(GO_DATAfolder)){
    dir.create(file.path(GO_DATAfolder))
  }
  exist_GOdata =  list.files(file.path(GO_DATAfolder))
  GO_DATAfile = paste(OrgDb_name,"_",OrgDb_version,"_",keyType,".RDS",sep = "")

  GO_DATA_environment = get_GOdata_Env()

  if(GO_DATAfile %in% exist_GOdata){
    use_cashed = F  ## load a exiting RDS
    if(exists("GO_DATAfile",envir = GO_DATA_environment) & exists("GO_DATA",envir = GO_DATA_environment)){
      GO_DATAfile2 = get("GO_DATAfile",envir = GO_DATA_environment)
      if(GO_DATAfile2 == GO_DATAfile)
        use_cashed = T  ## a correct RDS was loaded
      # GO_DATA = get("GO_DATA",envir = GO_DATA_environment)
    }

    if(use_cashed){
      GO_DATA = get("GO_DATA",envir = GO_DATA_environment)
    }else{
      GO_DATA = readRDS(paste(GO_DATAfolder,"/",GO_DATAfile,sep = ""))
      assign("GO_DATA",GO_DATA,envir = GO_DATA_environment)
      assign("GO_DATAfile",GO_DATAfile,envir = GO_DATA_environment)

      cat("loaded a existing GOdata\n")
    }

  }else{
    ## creat a GO_DATA for particular GO options: OrgDb,ont,keyType
    GO_DATA <- list(ALL = as.list(get_GO_data(OrgDb, "ALL", keyType)),
                    BP = as.list(get_GO_data(OrgDb, "BP", keyType)),
                    MF = as.list(get_GO_data(OrgDb, "MF", keyType)),
                    CC = as.list(get_GO_data(OrgDb, "CC", keyType)))
    ## save the GO_DATA as an rds file
    saveRDS(GO_DATA,file = paste(GO_DATAfolder,"/",GO_DATAfile,sep = ""))

    assign("GO_DATA",GO_DATA,envir = GO_DATA_environment)
    assign("GO_DATAfile",GO_DATAfile,envir = GO_DATA_environment)

    cat(paste("saved a new GOdata"))
  }

  GO_DATA = as.environment(GO_DATA[[ont]])
  # cat(ont)
  # cat("GO_DATA finished\n\n")
  # GO_DATA <- get_GO_data(OrgDb, ont, keyType)

  if (missing(universe))
    universe <- NULL

  if (ont == "ALL" && !pool) {
    lres <- lapply(c("BP", "CC", "MF"), function(ont)
      suppressMessages(enrichGO(gene, OrgDb, keyType, ont,
                                pvalueCutoff, pAdjustMethod, universe,
                                qvalueCutoff, minGSSize, maxGSSize
      ))
    )

    lres <- lres[!vapply(lres, is.null, logical(1))]
    if (length(lres) == 0)
      return(NULL)

    df <- do.call('rbind', lapply(lres, as.data.frame))
    geneSets <- lres[[1]]@geneSets
    if (length(lres) > 1) {
      for (i in 2:length(lres)) {
        geneSets <- append(geneSets, lres[[i]]@geneSets)
      }
    }
    res <- lres[[1]]
    res@result <- df
    res@geneSets <- geneSets
  } else {
    res <- enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             USER_DATA = GO_DATA
    )

    if (is.null(res))
      return(res)
  }

  res@keytype <- keyType
  res@organism <- get_organism(OrgDb)
  if(readable) {
    res <- setReadable(res, OrgDb)
  }
  res@ontology <- ont

  if (ont == "ALL") {
    res <- add_GO_Ontology(res, GO_DATA)
  }
  return(res)
}

environment(enrichGO) = asNamespace("clusterProfiler")
assignInNamespace("enrichGO", value = enrichGO, ns = "clusterProfiler")



############# PPI function #############
library(data.table)
library(visNetwork)
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
  scriptdir <- system.file("PPIdata",package = "DEP2")
  version <- "11.0"

  speciesID = switch (speciesname,
                      "human" = "9606",
                      "mouse" = "10090",
                      "Arabidopsis thaliana" = "3702",
                      "rat" = "10116",
                      "Saccharomyces cerevisiae" = "4932",
                      "zebrafish" = "7955"
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

  strEnv <- get_string_Env()
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
  proteins <<- proteinsDf2
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

  proteins <<- get_proteins(proteinsDf)
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
  if(nrow(naDf) > 0 & !quiet) cat(paste("Warning:  we couldn't map to STRING ", as.integer((nrow(naDf)/nrow(tempDf))*100), "% of your identifiers" , sep=""))
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
  nwplot <- visNetwork(nodes,links3,width = "170%",height = "500px") %>%
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
