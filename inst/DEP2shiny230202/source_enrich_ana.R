#### for gene annotation -------



## to handle the shiny tree input, is canceled in this version
# shinytreeInput_treatment <- function(Treeinput){
#
#   ord = stringr::str_count(Treeinput, "\\.")
#   Treeinput2 = Treeinput[order(ord)]
#   Treeinput2 = strsplit(Treeinput2,"\\.")
#   dropitem = c()
#   for(i in 1:(length(Treeinput2)-1)){
#     contain_sub <- Treeinput2[-(1:i)] %>% sapply(., function(x){
#       all(Treeinput2[[i]] %in% x)
#     }) %>% any()
#     if(contain_sub) dropitem = c(dropitem,i)
#   }
#   Treeinput3 = Treeinput2[-dropitem]
#   Treeinput3 = sapply(Treeinput3, function(x){x[length(x)]})
#   return(Treeinput3)
# }


# *ORA functions -----------------------------------------------------------
kegg_list <- function(db){
  url <- paste0("http://rest.kegg.jp/list/", db, collapse = "")
  clusterProfiler:::kegg_rest(url)
}


ko2name <- function(ko){
  p <- clusterProfiler:::kegg_list("pathway")
  ko2 <- gsub("^ko", "path:map", ko)
  ko.df <- data.frame(ko = ko, from = ko2)
  res <- merge(ko.df, p, by = "from", all.x = TRUE)
  res <- res[, c("ko", "to")]
  colnames(res) <- c("ko", "name")
  return(res)
}

get_GOdata_Env <- function(){
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

  if(GO_DATAfile %in% exist_GOdata){ ## exit a RDS file
    use_cashed = F  ## load a exiting RDS
    if(exists("GO_DATAfile",envir = GO_DATA_environment) & exists("GO_DATA",envir = GO_DATA_environment)){
      GO_DATAfile2 = get("GO_DATAfile",envir = GO_DATA_environment)
      if(GO_DATAfile2 == GO_DATAfile)
        use_cashed = T  ## a correct RDS was loaded
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
    cat("creat new db")
    ## creat a GO_DATA for particular GO options: OrgDb,ont,keyType
    GO_DATA <- list(ALL = as.list(my_get_GO_data(OrgDb, "ALL", keyType)),
                    BP = as.list(my_get_GO_data(OrgDb, "BP", keyType)),
                    MF = as.list(my_get_GO_data(OrgDb, "MF", keyType)),
                    CC = as.list(my_get_GO_data(OrgDb, "CC", keyType)))

    ## save the GO_DATA as an rds file
    saveRDS(GO_DATA,file = paste(GO_DATAfolder,"/",GO_DATAfile,sep = ""))

    assign("GO_DATA",GO_DATA,envir = GO_DATA_environment)
    assign("GO_DATAfile",GO_DATAfile,envir = GO_DATA_environment)

    cat(paste("saved a new GOdata:",paste(GO_DATAfolder,"/",GO_DATAfile,sep = "")))
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
if(rlang::is_installed("clusterProfiler")){
  environment(enrichGO) = asNamespace("clusterProfiler")
  assignInNamespace("enrichGO", value = enrichGO, ns = "clusterProfiler")
}

# library(org.Hs.eg.db)
# clusterProfiler:::get_GO_data(OrgDb = org.Hs.eg.db, ont = "ALL",keytype = "ENTREZID")
# OrgDb = org.Hs.eg.db; ont = "ALL"; keytype = "ENTREZID"
my_get_GO_data <- function (OrgDb, ont, keytype)
{
  GO_Env <- clusterProfiler:::get_GO_Env()
  use_cached <- FALSE
  ont2 <- NULL
  if (exists("ont", envir = GO_Env, inherits = FALSE))
    ont2 <- get("ont", envir = GO_Env)
  if (exists("organism", envir = GO_Env, inherits = FALSE) &&
      exists("keytype", envir = GO_Env, inherits = FALSE) &&
      !is.null(ont2)) {
    org <- get("organism", envir = GO_Env)
    kt <- get("keytype", envir = GO_Env)
    if (org == clusterProfiler:::get_organism(OrgDb) && keytype == kt && (ont ==
                                                        ont2 || ont2 == "ALL") && exists("goAnno", envir = GO_Env,
                                                                                         inherits = FALSE)) {
      use_cached <- TRUE
    }
  }
  if (use_cached) {
    goAnno <- get("goAnno", envir = GO_Env)
    if (!is.null(ont2) && ont2 != ont) {
      goAnno <- goAnno[goAnno$ONTOLOGYALL == ont, ]
    }
  }else {
    OrgDb <- GOSemSim::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (!keytype %in% kt) {
      stop("keytype is not supported...")
    }
    kk <- keys(OrgDb, keytype = keytype)
    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    if (ont != "ALL") {
      goterms <- goterms[goterms == ont] ## get the go terms ID of certain ont from GO.db
    }

    library(DBI)
    pkgname = OrgDb$packageName
    dbfile <- system.file("extdata",
                          paste(sub(".db$","",pkgname) ,".sqlite",sep=""), package=pkgname)
    if (!file.exists(dbfile))
      stop("DB file '", dbfile, "' not found")
    con <- DBI::dbConnect(RSQLite::SQLite(), dbfile)
    # dbListTables(con) ## check tables in orgDb

    genestbl = switch(keytype,
                      "ENTREZID" = con %>% tbl("genes") %>% rename(gene_id = gene_id) %>% select("_id","gene_id"),
                      "SYMBOL" = con %>% tbl("gene_info") %>% rename(gene_id = symbol) %>% select("_id","gene_id"),
                      "ENSEMBL" = con %>% tbl("ensembl") %>% rename(gene_id = ensembl_id) %>% select("_id","gene_id") ,
                      "UNIPROT" = con %>% tbl("uniprot") %>% rename(gene_id = uniprot_id) %>% select("_id","gene_id")
    ) ## extract the geneID table
    # genestbl = con %>% tbl("genes")
    gos = names(goterms)
    gotbl = con %>% tbl("go_all") %>% filter(go_id %in% gos) ## extract the go terms table
    goAnno = left_join(gotbl,genestbl,by="_id") %>% select("gene_id","go_id") %>%
      as.data.table()
    # go2gene = split(gotbl$gene_id,gotbl$go_id) ## to list
    # goAnno = gotbl2[,c("gene_id","go_id")]
    goAnno <- goAnno[with(goAnno, order(go_id)),] %>%  na.omit(cols="gene_id") %>%
      unique() %>% as.data.frame()
    dbDisconnect(con)



    # goAnno <- stack(go2gene)
    colnames(goAnno) <- c(keytype, "GOALL")
    # goAnno <- unique(goAnno[!is.na(goAnno[, 1]), ])
    goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]
    assign("goAnno", goAnno, envir = GO_Env)
    assign("keytype", keytype, envir = GO_Env)
    assign("ont", ont, envir = GO_Env)
    assign("organism", clusterProfiler:::get_organism(OrgDb), envir = GO_Env)
  }
  GO2GENE <- unique(goAnno[, c(2, 1)])
  GO_DATA <- clusterProfiler:::build_Anno(GO2GENE, clusterProfiler:::get_GO2TERM_table()) #大概3秒
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  if (!is.null(ont2) && ont2 == "ALL") {
    return(GO_DATA)
  }
  goOnt <- goOnt.df[, 2]
  names(goOnt) <- goOnt.df[, 1]
  assign("GO2ONT", goOnt, envir = GO_DATA)
  return(GO_DATA)
}



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


the1stname <- function(gene.name){
  names <- strsplit(as.character(gene.name),";")[[1]][1]
  return(names)
}




##go analysis function for DEP

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "Human", "mouse" or "rat", ...
#must be form reat <- goAnalysis(), because following use name reat
#df_with_lg2fc: if df_with_lg2fc
#species_df: data frame which have columns species, pkg and organism
goAnalysis <- function( gene_id, organism="Human", species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  # ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  # ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  #
  # ids <- inner_join(ids2, df, by = "name")
  #
  # if(df_with_lg2fc){
  #   ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  #   de = ids$fc
  #   names(de) = unlist(ids$ENTREZID)
  #   de = sort(de, decreasing = T)
  # } else {
  #   ids <- ids[!is.na(ids$ENTREZID), ]
  #   de = unlist(ids$ENTREZID)
  #   names(de) = de
  # }
  reat_ALL <- clusterProfiler::enrichGO(gene = names(gene_id), OrgDb = orgDB, ont = "ALL",
                                        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

  reat_BP <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "BP")
  reat_BP@ontology = "BP"

  reat_CC <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "CC")
  reat_CC@ontology = "CC"

  reat_MF <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "MF")
  reat_MF@ontology = "MF"

  reat <- list(ALL = reat_ALL, BP = reat_BP, CC = reat_CC, MF = reat_MF, de = gene_id)
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
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue <= pCutoff, p.adjust <= p.adj.cutoff, qvalue <= q.cutoff)

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
if(rlang::is_installed("enrichplot")){
  environment(Heatplot) = asNamespace("enrichplot")
}


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
if(rlang::is_installed("enrichplot")){
  environment(Cnetplot) = asNamespace("enrichplot")
  environment(fc_palette) = asNamespace("enrichplot")
}
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
if(rlang::is_installed("enrichplot")){
  environment(Emapplot) = asNamespace("enrichplot")
}

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
# library(clusterProfiler)
library(ggplot2)
# library(org.Mm.eg.db)
# library(org.Hs.eg.db)
# library(org.Rn.eg.db)
# library(DOSE)
# library(openxlsx)
# library(GO.db)
library(dplyr)
# library(ReactomePA)

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "Human", "mouse" or "rat", ...
#must be form reat <- keggAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
keggAnalysis <- function( gene_id ,
                          # df,
                          organism="Human", df_with_lg2fc = FALSE, species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  organism <- species_df$organism[species_df$species == organism]

  reat <- enrichKEGG(gene = names(gene_id), organism = organism,
                     pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = gene_id)
  return(reat)
}

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism:  one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
#must be form reat <- keggAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
# return a enrichResult
reactAnalysis <- function(gene_id, organism="Human", df_with_lg2fc = FALSE, species_df){
  pkg = species_df$pkg[species_df$species == organism]
  print(pkg)
  # cat("a0")
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  organism <- species_df$organism[species_df$species == organism]


  # cat("a1")
  reat <- ReactomePA::enrichPathway(gene = names(gene_id), organism = organism,
                                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable=TRUE)
  # cat("a4")
  reat <- list(res = reat, de = de)
  return(reat)
}

# reat: keggAnalysis or reactAnalysis funtion return value
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# q.cutoff: qvalue cutoff
give_res_and_table <- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.05, q.cutoff = 0.2){
  res <- reat$res
  sig_res <- clusterProfiler.dplyr::filter(res, pvalue <= pCutoff, p.adjust <= p.adj.cutoff, qvalue <= q.cutoff)

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
emaplot_for <- function(res, ShowCategory = 30, color = "p.adjust", layout = "kk") {
  try(Emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = layout), silent = T)
}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# circular: if circular layout
# colorEdge: if colorEdge
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
cnetplot_for <- function(res, ShowCategory = 5, circular = TRUE, colorEdge = TRUE, df_with_lg2fc = FALSE) {
  if(df_with_lg2fc) {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange = res$de, circular = circular, colorEdge = colorEdge), silent = T)
  } else {
    try(Cnetplot(x = res$sig_res, showCategory = ShowCategory, foldChange=NULL, circular = circular, colorEdge = colorEdge), silent = T)
  }

}

# res: givekegg_reat_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# df_with_lg2fc: if df_with_lg2fc , the same as function goAnalysis
heatplot_for <- function(res, ShowCategory = 30, df_with_lg2fc = FALSE) {
  if(df_with_lg2fc) {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = res$de), silent = T)
  } else {
    try(Heatplot(res$sig_res, showCategory = ShowCategory, foldChange = NULL), silent = T)
  }
}



# *GSEA functions ----------------------------------------------------------

########## gsego analysis function for DEP
#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "Human", "mouse" or "rat", ...
#must be form reat <- goAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
gsegoAnalysis <- function(
  # df,
  gene_list,
  organism="Human", species_df,
  topN = NULL){

  # cat(paste0("theOrg:", organism))
  pkg = species_df$pkg[species_df$species == organism]
  # cat(paste0("thePKG:", pkg))
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  if(!is.null(topN)){
    gene_list = gene_list[order(abs(gene_list),decreasing = T)[1:min(topN,length(gene_list))]]
  }
  # ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  # ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  #
  # ids <- inner_join(ids2, df, by = "name")
  #
  # ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  # de = ids$fc
  # names(de) = unlist(ids$ENTREZID)
  # de = sort(de, decreasing = T)

  set.seed(10086)

  reat_ALL <- suppressWarnings(try(gseGO(gene = gene_list, OrgDb = orgDB, ont = "ALL",
                                         pAdjustMethod = "BH",
                                         nPerm = 1000, ## use fgsea
                                         pvalueCutoff = 1, minGSSize = 3, verbose = FALSE, seed = FALSE), silent = TRUE))
  reat_ALL <- setReadable(reat_ALL, OrgDb = orgDB, keyType="ENTREZID")

  reat_BP <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "BP")
  reat_BP@setType = "BP"

  reat_CC <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "CC")
  reat_CC@setType = "CC"

  reat_MF <- clusterProfiler.dplyr::filter(reat_ALL, ONTOLOGY == "MF")
  reat_MF@setType = "MF"

  reat <- list(ALL = reat_ALL, BP = reat_BP, CC = reat_CC, MF = reat_MF, de = gene_list)
  return(reat)
}

# reat: gsegoAnalysis funtion return value
# ont: one of ALL, BP, CC, MF
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# NES.cutoff: NES cutoff,it means abs(NES) > NES.cutoff
# simplify: if remove redundancy of enriched GO terms
# Phenotype: the Phenotype you want to show, one or both of c("activated", "suppressed"), default: c("activated", "suppressed")
give_gseGO_res_and_table <- function(reat, ont = "BP", pCutoff = 0.05, p.adj.cutoff = 0.25, NES.cutoff = 1, simplify = FALSE, Phenotype = c("activated", "suppressed")){
  res <- reat[[ont]]
  res@result$phenotype = ifelse(res@result$NES > 0, "activated", "suppressed")
  res@result$group = paste(res@result$ONTOLOGY, res@result$phenotype, sep = "_")

  if(length(Phenotype) == 1) {
    res <- clusterProfiler.dplyr::filter(res, phenotype == Phenotype)#
  }

  sig_res <- clusterProfiler.dplyr::filter(res, pvalue <= pCutoff, p.adjust <= p.adj.cutoff, abs(NES) >= NES.cutoff)

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

# res: give_gseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
gse_dotplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(Split) {
      enrichplot:::dotplot(res, showCategory = ShowCategory, color = color, x = "NES", split = "group") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free")
    } else {
      enrichplot:::dotplot(res, showCategory = ShowCategory, color = color, x = "NES", split = "phenotype") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free")
    }
  } else {
    enrichplot:::dotplot(res, showCategory = ShowCategory, color = color, x = "NES", split="group")
  }
}

# res: give_gseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# Split: act when ont == "ALL", if split by ontology , TRUE or FALSE
gse_barplot <- function(res, ShowCategory = 20, color = "p.adjust", ont = "BP", Split = FALSE) {
  library(clusterProfiler.dplyr)
  library(ggplot2)
  if(ont == "ALL") {
    if(Split) {
      enrichplot:::barplot.enrichResult(res, showCategory = ShowCategory, color = color, x = "NES", split = "group") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free") + labs(y = "NES")
    } else {
      enrichplot:::barplot.enrichResult(res, showCategory = ShowCategory, color = color, x = "NES", split = "phenotype") + facet_grid(ONTOLOGY ~ ., scale="free", space = "free") + labs(y = "NES")
    }
  } else {
    # barplot.enrichResult
    enrichplot:::barplot.enrichResult(res, showCategory = ShowCategory, color = color, x = "NES", split="group") + labs(y = "NES")
  }
}

# res: giveGO_res_and_table or give_gseGO_res_and_table funtion return value
# ShowCategory: number of go terms of plot
# color: p.adjust" or "pvalue"
# try(emapplot(res$sig_res, showCategory = ShowCategory, color = color, layout = "kk"), silent = T)


# res: giveGO_res_and_table or give_gseGO_res_and_table funtion return value
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


#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "Human", "mouse" or "rat", ...
#must be form reat <- gsekeggAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
gsekeggAnalysis <- function(gene_list, organism="Human", species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  organism <- species_df$organism[species_df$species == organism]

  set.seed(1234)
  reat <- suppressWarnings(try(gseKEGG(gene = gene_list, organism = organism,
                                       pAdjustMethod = "BH",
                                       nPerm = 1000,
                                       pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE))
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = de)
  return(reat)
}

#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism:  one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
#must be form reat <- gsereactAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
gsereactAnalysis <- function(gene_list, organism="Human", species_df){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  organism <- species_df$reactome_organism[species_df$species == organism]

  # ids1 = my_to_entrezid(orgDB = orgDB, gene = as.character(df$name))
  # ids2 <- ids1 %>% tibble::rownames_to_column() %>% dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)
  #
  # ids <- inner_join(ids2, df, by = "name")
  #
  # ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]
  # de = ids$fc
  # names(de) = unlist(ids$ENTREZID)
  # de = sort(de, decreasing = T)

  set.seed(1234)
  reat <- suppressWarnings(try(ReactomePA::gsePathway(gene = gene_list, organism = organism,
                                                      pAdjustMethod = "BH",
                                                      nPerm = 1000,
                                                      pvalueCutoff = 1, verbose = FALSE, seed = FALSE), silent = TRUE))
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = de)
  return(reat)
}


# reat: gsekeggAnalysis or gsereactAnalysis funtion return value
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# NES.cutoff: NES cutoff,it means abs(NES) > NES.cutoff
# Phenotype: the Phenotype you want to show, one or both of c("activated", "suppressed"), default: c("activated", "suppressed")
give_gsegsereat_res_and_table <- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.25, NES.cutoff = 1, Phenotype = c("activated", "suppressed")){
  res <- reat$res
  res@result$phenotype = ifelse(res@result$NES > 0, "activated", "suppressed")

  if(length(Phenotype) == 1) {
    res <- clusterProfiler.dplyr::filter(res, phenotype == Phenotype)#
  }

  sig_res <- clusterProfiler.dplyr::filter(res, pvalue <= pCutoff, p.adjust <= p.adj.cutoff, abs(NES) >= NES.cutoff)

  res_table <- as.data.frame(res)
  sig_res_table <- as.data.frame(sig_res)

  table <- list(all_table = res_table, sig_table = sig_res_table, sig_res = sig_res, de = reat[["de"]])

  return(table)
}

# reat: gsekeggAnalysis or gsereactAnalysis funtion return value
# pCutoff: pvalue cutoff
# p.adj.cutoff: padj cutoff
# NES.cutoff: NES cutoff,it means abs(NES) > NES.cutoff
# Phenotype: the Phenotype you want to show, one or both of c("activated", "suppressed"), default: c("activated", "suppressed")
# give_gsekegg_gsereat_res_and_table <-
give_gseOther_res_and_table<- function(reat, pCutoff = 0.05, p.adj.cutoff = 0.25, NES.cutoff = 1, Phenotype = c("activated", "suppressed")){
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





## $add for msigdb analysis


#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "Human", "mouse" or "rat", ...
#must be form reat <- keggAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
msigdb_Analysis <- function(gene_id, organism="Human", df_with_lg2fc = FALSE, species_df, category, subcategory){
  library(msigdbr)
  library(enrichplot)
  library(clusterProfiler)
  speciesID = switch (organism,
                      "Human" = "Homo sapiens",
                      "Mouse" = "Mus musculus",
                      "Rat" = "Rattus norvegicus")

  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  # organism <- species_df$organism[species_df$species == organism]
  #
  if(!is.null(category)){
    m_t2g = msigdbr(species = speciesID, category = category, subcategory = subcategory) %>% dplyr::select(gs_name, entrez_gene)
  }else{
    m_t2g = msigdbr(species = speciesID, category = NULL, subcategory = NULL) %>%
      filter(gs_cat %in% subcategory | gs_subcat %in% subcategory) %>% dplyr::select(gs_name, entrez_gene)
    m_t2g = m_t2g
  }
  reat <- enricher(gene = names(gene_id), TERM2GENE = m_t2g, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)

  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = gene_id)
  return(reat)
}





#df: left gene right fc or log2fc (colnames is always fc), the output legend is awlays fold change
#organism: "Human", "mouse" or "rat", ...
#must be form reat <- gsekeggAnalysis(), because following use name reat
#species_df: data frame which have columns species, pkg and organism
gsemsigdb_Analysis <- function(gene_list, organism="Human", species_df, category, subcategory
                               # msigdb = c("c2.cgp (chemical and genetic perturbations)", "c2.cp.biocarta", "c2.cp (Canonical pathways)", "c3.all (motif gene sets)", "c3.tft (transcription factor targets)", "c6.all (oncogenic signatures)", "c7.all (immunologic signatures)", "h.all.v6.1 (hallmark gene sets)")){

){
  speciesID = switch (organism,
                      "Human" = "Homo sapiens",
                      "Mouse" = "Mus musculus",
                      "Rat" = "Rattus norvegicus")

  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  if(!is.null(category)){
    m_t2g = msigdbr(species = speciesID, category = category, subcategory = subcategory) %>% dplyr::select(gs_name, entrez_gene)
  }else{
    m_t2g = msigdbr(species = speciesID, category = NULL, subcategory = NULL) %>%
      filter(gs_cat %in% subcategory | gs_subcat %in% subcategory) %>% dplyr::select(gs_name, entrez_gene)
    m_t2g = m_t2g
  }


  set.seed(1234)
  reat <- try(GSEA(geneList = gene_list, TERM2GENE = m_t2g, pAdjustMethod = "BH", pvalueCutoff = 1, verbose = TRUE, seed = FALSE), silent = TRUE)
  reat <- setReadable(reat, OrgDb = orgDB, keyType="ENTREZID")

  reat <- list(res = reat, de = de)
  return(reat)
}




