############# PPI function #############
# library(data.table)
# library(visNetwork)
# library(igraph)
# library(htmlwidgets)

## creat and call a environment .string_Env
get_string_Env <- function () {
  if (!exists(".string_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".string_Env", new.env(), envir=envir)
  }
  get(".string_Env", envir = .GlobalEnv)
}


# stringNetwork
## PPI analysis on candidates proteins/genes according the
#' Protein-protein interaction analysis
#'
#' \code{test_PPI} search protein-protein interaction on significant candidate
#' based on a local STRING database.
#'
#' @param x A SummarizedExperiment/DEGdata output from \code{\link{add_rejections}} or
#' a charachter vector containing candidate identifiers(SYMBOL, EntrezID, UniprotID or ENSEMBL).
#' @param contrasts Character, specifies the contrasts to get significant candidates.
#' @param species The species name.
#' @param STRING.version Character(1),the version of STRING data.
#' @param choose_scores NULL or "combined_score" or
#' a subset of c("neighborhood","fusion","cooccurence","coexpression",
#' "experimental","database","textmining"). If is NULL or "combined_score",
#' cutoff is on the scores combining above all evidence, otherwise is on the subset.
#' @param score_cutoff A numeric lager than 400, required lowest interaction scores.
#' A score lowwer than 400 means the interaction is unconfident
#'
#' @details
#' Identifiers in x are firstly transformed to entrezID and then mapped the relative STRING id.
#' Next the protein-protein interaction is exctracted from STRING.link.detail.
#' The local STRING data (including protein.info, protein.aliases and protein.links.detailed) is necessary for this function.
#' If it is absent, \code{test_PPI} will download from \href{http://stringdb-static.org/download/}{STRING} automatically according
#' \code{species} and \code{STRING.version}.
#'
#' @return
#' A data.frame of PPI information
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example
#' data(Silicosis_pg)
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Differential test
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols,mode = "delim")
#' filt <- filter_se(se, thr = 0, fraction = 0.4, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.05)
#' diff <- test_diff(imputed, type = "control", control  = c("PBS"), fdr.type = "Storey's qvalue")
#' dep <- add_rejections(diff, alpha = 0.01,lfc = 2)
#'
#' # PPI construct
#' if(check_PPI_depends()){
#'
#'   # Load STRING data. If local STRING is missing, PPI_res will download to system.file("PPIdata", "Mouse",package = "DEP2")
#'   load_PPIdata(speciesname = "Mouse")
#'
#'   PPI_res <- test_PPI(dep, contrasts = "W4_vs_PBS", species = "Mouse",
#'                       STRING.version = "11.5",score_cutoff = 400)
#'   head(PPI_res)
#' }
#' }
test_PPI <- function(x,
                     contrasts = NULL,
                     species = "Human",
                     STRING.version = "11.5",
                     choose_scores = NULL,
                     score_cutoff = 400
                     ){
  check_pak <- check_PPI_depends()
  if(!isTRUE(check_pak)){
    stop("Packages ",paste0(check_pak,collapse = ", ")," are required for PPI function, but not found")
  }

  assertthat::assert_that(class(x) == "SummarizedExperiment"|class(x) == "DEGdata"|class(x) == "character",
                          is.null(contrasts)| is.character(contrasts),
                          is.character(species) && length(species) ==1,
                          is.character(STRING.version) && length(STRING.version) == 1,
                          is.null(choose_scores)|
                            (is.character(choose_scores) &
                               all(choose_scores %in%
                                     c("combined_score","neighborhood","fusion","cooccurence","coexpression",
                                       "experimental","database","textmining"))
                            ),
                          is.numeric(score_cutoff)|is.integer(score_cutoff) )
  the_annoSpecies_df = annoSpecies_df()
  speciesID = the_annoSpecies_df$speciesID[which(the_annoSpecies_df$species == species)]

  if(class(x) == "character"){
    gene = x
  }else if(class(x) == "DEGdata" | class(x) == "SummarizedExperiment"){
    sig = get_signicant(x, contrasts = contrasts,return_type = "subset")
    if( all(c("gene_name","protein_ID") %in% colnames(rowData(x))) ){
      ## check may the object is a PTM analyses generated from make_unique_ptm
      gene = rowData(x)$gene_name ## use the gene name but not the PTM identifiers.
    }else{
      gene = rownames(sig) ## extract rownames
    }
  }

  if(length(gene) < 2){
    stop("Giving candidate is too few, please check your input")
  }

  orgDB = the_annoSpecies_df$pkg[the_annoSpecies_df$species == species]
  require(orgDB, character.only = TRUE)
  orgDB = get(orgDB)
  gene_id_table <- map_to_entrezid(gene, orgDB = orgDB)
  ids <- gene_id_table %>% tibble::rownames_to_column() %>%
    dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID) %>%
    filter(!is.na(ENTREZID)) %>% filter(!duplicated(ENTREZID))
  if(nrow(ids) < 2){
    stop("No giving genes/proteins map in the database, please check your input!")
    # return("No mapped genes")
  }


  message("Check and load STRING data...")
  message(load_PPIdata(speciesname = species, STRING.version = STRING.version))
  strEnv = get_string_Env()
  aliasDf <- get("aliasDf",envir = strEnv)
  proteinsDf <- get("proteinsDf",envir = strEnv)
  protein_links_detail <- get("protein_links_detail",envir = strEnv)


  message("Mapping STRING ID...")
  if(is.null(choose_scores) || choose_scores == ""){
    choose_scores = "combined_score"
  }else if(length(choose_scores)%in% 2:(ncol(protein_links_detail)-3) && "combined_score" %in% choose_scores){
    warning("choose_scores should be NULL or 'combined_score' or a subset of sub scores in c('neighborhood','fusion','cooccurence'",
            "',coexpression','experimental','database','textmining'). You have chosen 'combined_score' and the other choosen sub scores is invalid,",
            "all scores are used for combined_score")
    choose_scores = colnames(protein_links_detail)[-c(1,2)]
  }else if("combined_score" %in% choose_scores){
    choose_scores = colnames(protein_links_detail)[-c(1,2)]
  }
  data_mapped <- ids %>% mapPPI(my_data_frame = .,
                                   my_data_frame_id_col_names = "ENTREZID",
                                   aliasDf = aliasDf,
                                   proteinsDf = proteinsDf, ## data_mapped 3cols ENTREZID SYMBOL(name) STRING_id
                                   removeUnmappedRows = TRUE,reload = F)
  message("Getting interactions...")
  link_table = get_interactions(protein_detail = protein_links_detail,
                                mapped_data = data_mapped,choose_scores = choose_scores, score_cutoff = score_cutoff)

  links = link_table %>%                                                ## links 3cols from to weight from&to is SYMBOL
    mutate(from = data_mapped[match(from, data_mapped$STRING_id), "name"]) %>%
    mutate(to = data_mapped[match(to, data_mapped$STRING_id), "name"])

  ## filter out duplicate links and low scorecutoff
  links2 = filterlink(links, scorecutoff = score_cutoff)
  message("Finished")
  return(links2)
}

############# network function for PPI #############
##
#' Draw a network on PPI result
#'
#' Draw a network on PPI result output from \code{test_PPI}
#'
#' @param PPIlinks The PPI result from {test_ppi}
#' @param layoutway Character(1), the name of a layout function in
#' \code{igraph} package,such as "layout_as_star", "layout_as_tree",
#' "layout_in_circle". See \link[igraph]{layout}.
#' @param nodecolor Character(1), the color of node.
#' @param nodeshape Character(1), the shape of node. One of "circle", "square", "raster",
#' or other vertex.shapes in igraph.
#' @param linecolor Character(1), the color of edge line.
#' @param nodesize Numeric(1), the size of node.
#' @param changesize Logical(1), whether change node size by link number.
#' @param fontsize Numeric(1), the font size of label.
#' @param linewidth Numeric(1), the width of edge.
#' @param changewidth Logical(1), whether change edge width by interaction score.
#' @param smoothline Logical(1), whether smooth the edge. Works only for visNetwork
#' @param smoothtype Character(1), one of 'dynamic', 'continuous', 'discrete', 'diagonalCross',
#' 'straightCross', 'horizontal', 'vertical', 'curvedCW', 'curvedCCW', 'cubicBezier'
#' @param highlightkey logical(1), whether label the central nodes by red
#' which have at least 6 interactions.
#' @param returntype Character(1), "visNetwork" or "igraph".
#'
#' @return
#' A visNetwork plot or a igraph obejct according \code{returntype}.
#'
#' @examples
#' \dontrun{
#' # Load example
#' data(Silicosis_pg)
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Differential test
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols,mode = "delim")
#' filt <- filter_se(se, thr = 0, fraction = 0.4, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.05)
#' diff <- test_diff(imputed, type = "control", control  = c("PBS"), fdr.type = "Storey's qvalue")
#' dep <- add_rejections(diff, alpha = 0.01,lfc = 2)
#'
#' # PPI construct
#' if(check_PPI_depends()){
#'
#'   # Load STRING data. If local STRING is missing, PPI_res will download to system.file("PPIdata", "Mouse",package = "DEP2")
#'   load_PPIdata(speciesname = "Mouse")
#'
#'   PPI_res <- test_PPI(dep, contrasts = "W4_vs_PBS", species = "Mouse",
#'                       STRING.version = "11.5",score_cutoff = 400)
#'
#'   ## igraph network
#'   PPI_ig <- PPInetwork(PPI_res, returntype = "igraph")
#'   igraph::plot.igraph(PPI_ig)
#'
#'   ## visNetwork network
#'   PPInetwork(PPI_res, returntype = "visNetwork")
#'   }
#' }
#'
#' @export
PPInetwork <- function(PPIlinks, layoutway = "layout_components",nodecolor = "#2EA9DF",nodeshape = c("circle","square"),linecolor = "#ADD8E6",
                          nodesize = 20,changesize=FALSE,
                          fontsize=25,changewidth=FALSE,linewidth = 5, smoothline=FALSE,smoothtype="continous",highlightkey=TRUE,
                          returntype = c("visNetwork","igraph")){
  # PPIlinks = linksTable
  nodeshape = match.arg(nodeshape)
  returntype = match.arg(returntype)
  links3 <- PPIlinks
  nodes <- PPIlinks %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
  net <- igraph::graph_from_data_frame(d=links3,vertices=nodes,directed = F)
  igraph::graph_attr(net, "layout") <- get(layoutway, envir = asNamespace("igraph"))
  if(nodeshape == "circle"){
    nodes$shape = "dot"
  }else{
    nodes$shape = nodeshape
  }

  igraph::V(net)$shape <- nodeshape

  if(changewidth){
    links3$width <- links3$combined_score/1500*linewidth
    igraph::E(net)$width = links3$combined_score/1500*linewidth
  }else{
    links3$width =linewidth
    igraph::E(net)$width = linewidth/4
  }

  if(changesize){
    nodes$size <- igraph::degree(net)*nodesize/2
    igraph::V(net)$size <- igraph::degree(net)*nodesize/2/6
    # cat(changesize)
  }else{
    nodes$size = nodesize
    igraph::V(net)$size <- nodesize/6
  }
  if(highlightkey){
    nodes$color <- ifelse(igraph::degree(net)>=6,"#B54434",nodecolor)
    igraph::V(net)$color <- ifelse(igraph::degree(net)>=6,"#B54434",nodecolor)
  }else{
    nodes$color = nodecolor
    igraph::V(net)$color <- nodecolor
  }
  igraph::V(net)$label.cex <- fontsize/50
  igraph::V(net)$label.dist <- -0.17*igraph::V(net)$size

  # nodes$size =200
  colnames(nodes)[1] <- "id"

  nwplot <- visNetwork::visNetwork(nodes,links3,width = "170%",height = "500px") %>%
    visNetwork::visIgraphLayout(layout = layoutway ) %>%
    visNetwork::visNodes(size=nodes$size,
             # color = nodecolor ,
             font= list(size= fontsize,align = "left"),
             shadow = list(enabled = TRUE, size = 10)) %>%
    visNetwork::visEdges(shadow = F,
             color = list(color = linecolor, highlight = "red"),
             smooth = list(enabled=smoothline,type=smoothtype)
    ) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = T, hover = T),
               nodesIdSelection = T)%>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
               forceAtlas2Based = list(gravitationalConstant = -500))

  if(returntype == "visNetwork"){
    return(nwplot)
  }else if(returntype == "igraph"){
    return(net)
  }
}


# download a file when it is not already present. From STRINGdb package.
downloadAbsentFile <- function(urlStr, oD = tempdir()){
  fileName = tail(strsplit(urlStr, "/")[[1]], 1)
  temp <- paste(oD,"/", fileName,sep="")
  temp2 <- paste0(temp,".tmp")
  fristDownload = F
  if(! file.exists(temp) || file.info(temp)$size==0){
    cat("Download ", fileName," to ",oD,"\n")
    # download.file(urlStr,temp2, timeout = 400)
    downloader::download(urlStr,
                         temp2, timeout = 400)
    # curl::curl_download(urlStr,
    #                     temp2,quiet = FALSE)
    file.rename(temp2,temp)
    fristDownload = T
  }
  if(file.info(temp)$size==0) {
    unlink(temp)
    temp=NULL
    cat(paste("ERROR: failed to download ", fileName,".\nPlease check your internet connection and/or try again. " ,sep=""))
  }
  cat("fristDownload\n")
  return(fristDownload)
}



Check_STRINGdata <- function(speciesID, STRING.version = "11.5", STRINGdata_path = NULL){
  file_version = as.character(STRING.version)
  protocol <- "https"
  proteinsDffilename = paste(protocol,
                             "://stringdb-static.org/download/protein.info.v",
                             file_version, "/", speciesID, ".protein.info.v", file_version, ".txt.gz", sep="")
    # paste(speciesID, ".protein.info.v", version, ".txt.gz", sep="")
  aliasesfilename = paste(protocol,
                          "://stringdb-static.org/download/protein.aliases.v",
                          file_version, "/", speciesID, ".protein.aliases.v", file_version, ".txt.gz", sep="")
    # paste(speciesID, ".protein.aliases.v", version, ".txt.gz", sep="")
  links_detailfilename = paste(protocol,
                               "://stringdb-static.org/download/protein.links.detailed.v",
                               file_version, "/", speciesID, ".protein.links.detailed.v", file_version, ".txt.gz", sep="")
    # paste(speciesID, ".protein.links.detailed.v", version, ".txt.gz", sep="")

  if(is.null(STRINGdata_path)){
    STRINGdata_path <-  file.path(system.file(package = "DEP2"),"PPIdata",speciesname)
  }

  if(!dir.exists(STRINGdata_path)){
    dir.create(STRINGdata_path, recursive = T)
    cat(paste("PPIdata directory is",STRINGdata_path,"\n"))
  }
  A = downloadAbsentFile(proteinsDffilename,oD = STRINGdata_path)
  B = downloadAbsentFile(aliasesfilename,oD = STRINGdata_path)
  C = downloadAbsentFile(links_detailfilename,oD = STRINGdata_path)

  if(C){
    cat("STRING data version",file_version, "for species", speciesID, "downloaded finished\n")
    trimSTRINGdata(links_detailfilename, oD = STRINGdata_path)
  }
}

## trim the firstly download link.detail file and replace it to save storage space
trimSTRINGdata <- function(links_detailfilename, oD = NULL)  {
  links_detailfilename = tail(strsplit(links_detailfilename, "/")[[1]], 1)
  links_detailfile = file.path(oD, links_detailfilename)
  if(!file.exists(links_detailfile))
    stop("links.detailed file not exists")
  links_detail = fread(links_detailfile)
  ## just save links_detailed which combined_score >= 400
  links_detail = links_detail %>% filter(combined_score >= 400)
  message("Trime ", links_detailfilename, " by combined score > 400.")
  data.table::fwrite(links_detail, file = links_detailfile)
}


## load PPIdata of spcies from DEP/PPIdata dir to strEnv=.string_Env
#'
#' Check and load a local STRING data. If local file do not exist, will try to download from \href{http://stringdb-static.org/download/}{STRING}.
#' @param speciesname Species name.
#' @param STRING.version character(1), STRING dataset version, default is "11.5".
#' @param STRINGdata_path NULL or a selected path. Recommend NULL.
#'
#' @importFrom data.table fread
#' @export
#' @examples
#' \dontrun{
#' # Load STRING data. If local STRING is missing, PPI_res will download to system.file("PPIdata", "Mouse",package = "DEP2")
#' load_PPIdata(speciesname = "Mouse")
#' }
load_PPIdata <- function(speciesname, STRING.version = "11.5", STRINGdata_path = NULL){
  strEnv <- get_string_Env()

  the_annoSpecies_df <- annoSpecies_df()
  speciesID = the_annoSpecies_df$speciesID[which(the_annoSpecies_df$species == speciesname)]
  STRING.version <- as.character(STRING.version)

  # check whether STRING data was loaded to strEnv
  if(exists(c("load succed"),envir = strEnv) && get("speciesname",envir = strEnv) == speciesname &&
     get("load succed",envir = strEnv) && get("STRING.version",envir = strEnv) == STRING.version){
    return("use_cache")
  }

  if(is.null(STRINGdata_path)){
    STRINGdata_path <- file.path(system.file(package = "DEP2"),"PPIdata",speciesname)
  }

  cat(paste("Input species is ",speciesname ,"\t", "The species ID is ",speciesID,"\n",sep=""))
  Check_STRINGdata(speciesID, STRING.version, STRINGdata_path = STRINGdata_path)

  ### Load STRING files to string_Env
  proteinsDffile <- paste(STRINGdata_path, "/", speciesID, ".protein.info.v", STRING.version, ".txt.gz", sep="")
  aliasesfile = paste(STRINGdata_path, "/", speciesID, ".protein.aliases.v", STRING.version, ".txt.gz", sep="")
  links_detailfile = paste(STRINGdata_path, "/", speciesID, ".protein.links.detailed.v", STRING.version, ".txt.gz", sep="")


  if(!all( file.exists(c(proteinsDffile,aliasesfile,links_detailfile)) )){
    stop(paste("Lack the necessary String datasets in ",file.path(STRINGdata_path,speciesname),
               "\n Please download the info,aliases,detailed txt.gz files of this species from String, and paste to ",
               file.path(STRINGdata_path,speciesname),sep = ""))
  }

  message("load STRING data from ",paste0(STRINGdata_path,"/"))
  message("lood protein.info")
  assign("proteinsDf",data.table::fread(proteinsDffile, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE, quote=""),envir = strEnv)
  message("lood protein.aliases")
  assign("aliasDf",suppressWarnings(data.table::fread(aliasesfile)),envir = strEnv)
  message("lood protein.links.datailed")
  assign("protein_links_detail",data.table::fread(links_detailfile),envir = strEnv)
  assign("speciesID",speciesID,envir= strEnv)
  assign("speciesname",speciesname,envir= strEnv)
  assign("STRING.version",STRING.version,envir= strEnv)
  assign("load succed",T,envir=strEnv)
  cat("load succeed")
  return("Loading finished!")
}


## get the interaction score table based on mapped_data.
get_interactions = function(protein_detail,mapped_data,choose_scores,score_cutoff = 400){
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


## modified from STRINGdb
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
  colnames(proteinsDf)[1] = "protein_external_id"
  proteinsDf2 = subset(proteinsDf, select=c("protein_external_id",  "preferred_name", "protein_size", "annotation"))
  proteins <- proteinsDf2
  # cat(paste("load proteinsDFfile: ",temp ,"\n", sep=""))

  return(proteins)
}

## map aliases of genes, modified from STRINGdb
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

# mapping function (add the possibility to map using more than one column), modified from STRINGdb::multi_map_df
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


## Map a STRING identifiers, modified from STRINGdb
mapPPI = function(my_data_frame,
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
  if(nrow(naDf) > 0 & !quiet) cat(paste("We couldn't map to STRING: ", as.integer((nrow(naDf)/nrow(tempDf))*100), "% of your identifiers" , sep=""))
  if(removeUnmappedRows) tempDf = subset(tempDf, !is.na(STRING_id))

  return(tempDf)
}




filterlink <- function(linkstable,scorecutoff = NULL){
  if(!is.null(scorecutoff))
    linkstable <- linkstable %>% dplyr::filter(combined_score>=scorecutoff)
  # drop = c()
  linkslist = linkstable[,c("from","to")] %>% apply(., 1, function(x){unname(sort(x))},simplify = F) ## remove the replicates like (from A to B) and (from B to A)
  drop = linkslist %>% duplicated()
  table(drop)
  linkstable = linkstable[!drop,]
  # for (i in 1:(nrow(linkstable)-1)) {
  #   if(i >= nrow(linkstable) ) break
  #   for (j in (i+1):nrow(linkstable)) {
  #     if(linkstable[i,1]==linkstable[j,2] &&linkstable[i,2]==linkstable[j,1]){
  #       linkstable = linkstable[-j,]
  #       break
  #     }
  #   }
  # }
  return(linkstable)
}

# linksTable = PPI_table; layoutway = "layout_components";nodecolor = "#2EA9DF";nodeshape = "dot";linecolor = "#ADD8E6";nodesize = 20;changesize=FALSE;
# fontsize=200;changewidth=FALSE;linewidth = 5;smoothline=FALSE;smoothtype="continous";highlightkey=TRUE


