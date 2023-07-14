`%<>%` <- magrittr::`%<>%`



#' Species information table
#'
#' Record the species names IDs and their org database package.
#'
#' @return
#' a data.frame
#'
#' @export
#'
#' @examples
annoSpecies_df <- function(){
  annoSpecies_df <-
    data.frame(
      species = c(
        "", "Anopheles", "Arabidopsis", "Bovine", "Worm",
        "Canine", "Fly", "Zebrafish", "E coli strain K12",
        "E coli strain Sakai", "Chicken", "Human", "Mouse",
        "Rhesus", "Malaria", "Chimp", "Rat",
        "Yeast", "Streptomyces coelicolor", "Pig", "Toxoplasma gondii",
        "Xenopus"
      ),
      pkg = c(
        "", "org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
        "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
        "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
        "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
        "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db",
        "org.Xl.eg.db"
      ),
      stringsAsFactors = FALSE
    )
  annoSpecies_df$organism = c("", "aga", "ath", "bta", "cel", "cfa", "dme", "dre", "eco", "ecs", "gga", "hsa", "mmu", "mcc", "pfa", "ptr", "rno", "sce", "sco", "ssc", "tgo", "xla")

  annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species), ]
  annoSpecies_df$reactome_organism = c("","anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human",
                                       "malaria", "mouse", "pig", "rat", "rhesus", "coelicolor", "gondii", "celegans", "xenopus", "yeast", "zebrafish")
  annoSpecies_df$msigdbr_species = c("","","","Bos taurus", "Canis lupus familiaris", "Gallus gallus", "Pan troglodytes",
                                     "","","Drosophila melanogaster", "Homo sapiens","", "Mus musculus", "Sus scrofa",
                                     "Rattus norvegicus", "Macaca mulatta", "", "", "Caenorhabditis elegans", "Xenopus tropicalis",
                                     "Saccharomyces cerevisiae", "Danio rerio")

  string_spe <- c( "Anopheles" = "7165", "Arabidopsis" = "3702", "Bovine" = "9913",
                   "Worm" = "6239", "Canine" = "9612", "Fly" = "7227", "Zebrafish" = "7955",
                   "E coli strain K12" = "511145", "Chicken" = "9031", "Human" = "9606",
                   "Mouse" = "10090", "Rhesus" = "9544", "Malaria" = "5833", "Chimp" = "9598",
                   "Rat" = "10116", "Yeast" = "4932", "Streptomyces coelicolor" = "100226",
                   "Pig" = "9823", "Toxoplasma gondii" = "5811", "Xenopus" = "8364")
  annoSpecies_df$speciesID = string_spe[match(annoSpecies_df$species, names(string_spe))]
  annoSpecies_df$speciesID[is.na(annoSpecies_df$speciesID)] = ""
  rownames(annoSpecies_df) <- annoSpecies_df$species

  ## filter out some unusual species
  annoSpecies_df = annoSpecies_df[which(annoSpecies_df$msigdbr_species != ""),]
  species_order = c("Human","Mouse", "Rat", "Yeast", "Pig", "Zebrafish")
  species_order = c(species_order, annoSpecies_df$species[!annoSpecies_df$species %in% species_order])
  annoSpecies_df = annoSpecies_df[match(species_order,annoSpecies_df$species),]

  return(annoSpecies_df)
}

annoSpecies_df2 <- function(){
  theannoSpecies_df <- DEP2:::annoSpecies_df()
  pkg <- theannoSpecies_df$pkg
  # installed <- pkg[-1] %>% sapply(.,function(x){rlang::is_installed(x)})
  # theannoSpecies_df[1+which(installed),]
  installed <- pkg %>% sapply(.,function(x){rlang::is_installed(x)})
  theannoSpecies_df[which(installed),]
}
# DEP2:::annoSpecies_df()

#' ID transform SE or DEGdata
#'
#' Transform ID for SummarizedExperiment or DEGdata according origin ID.
#' The annotation package of certain species must be installed. Using \code{annoSpecies_df}
#' to check species names and the required packages.
#'
#' @param x SummarizedExperiment object from \code{\link{make_se}()}) or \code{\link{make_pe}()}),
#' or a DEGdata object from \code{\link{test_diff_deg}()}).
#' @param from_columns Character(), the origin ID from one of "rownames" or column of \code{rowData(x)}
#' @param fromtype Character(1), the type of origin ID, e.g. "ENSEMBEL", "SYMBOL", "UNIPROT", "ENTREZID".
#' @param species Character(1), the species name.
#' @param replace_rowname NULL or character. Should be one of NULL, "SYMBOL", "ENTREZID", "UNIPROT", "ENSEMBL".
#' If it is not NULL, replace rownames(x) by selected kind of ID.
#'
#' @return
#' A SummarizedExperiment or DEGdata which rowData includes gene information like
#' "ENSEMBEL", "SYMBOL", "UNIPROT" or "ENTREZID" transformed from giving ID.
#' If x is a DEGdata, gene information table is also store in the geneinfo slot.
#' Rownames of object are replaced according replace_rowname.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(Silicosis_pg)
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols,mode = "delim")
#'
#' ## transform uniport ID
#' check_organismDB_depends(organism = "Mouse") # check annotation package of Mouse
#' head(rowData(se)$ID)
#' se2 <- ID_transform(se,from_columns = "ID",fromtype = "UNIPROT",species = "Mouse")
#' head(rowData(se2)$SYMBOL)
#' head(rownames(se2))
#' }
ID_transform <- function(x, from_columns = "rownames", fromtype = "ENSEMBL", species = "Human", replace_rowname = "SYMBOL")
{
  ## check org db
  check_pak <- check_organismDB_depends(organism = species, install = T)
  if(!isTRUE(check_pak)){
    stop("Packages ",paste0(check_pak,collapse = ", ")," are required for ID transform but not found")
  }

  if(from_columns == "rownames"){
    ids = rownames(x)
  }else{
    rd = rowData(x)
    ids = rd[, from_columns]
  }

  Species_df = DEP2:::annoSpecies_df()

  anno_db = Species_df[species, ]$pkg
  cat(anno_db,"\n")
  require(anno_db, character.only = TRUE)
  anno_db = get(anno_db)
  # anno_db = eval(parse(text = paste0(anno_db,"::",anno_db)))

  columns <- keytypes(anno_db)
  columns <- intersect(c("SYMBOL", "ENTREZID", "UNIPROT", "ENSEMBL"), columns)
  # columns <- columns[which(columns != input$idtype)]
  ann <- try(AnnotationDbi::select(anno_db, keys = ids, column = columns, keytype = fromtype, multiVals = "first"))
  if(class(ann) == "try-error"){
    stop("Mapping ID failed, please check your input")
  }

  ann = ann[!duplicated(ann[, fromtype]), ]
  ann2 = ann[match(ids,ann[, fromtype]),]

  rownames(x) = ifelse(!is.na(ann2[,replace_rowname]),make.unique(ann2[,replace_rowname]),rownames(x))

  if(class(x) == "DEGdata"){
    x@geneinfo = ann2
    if( !is.null(x@test_result) && nrow(x@test_result) == nrow(x)){
      x@test_result$symbol = ann2$SYMBOL
      x@test_result =  x@test_result %>% data.frame() %>%
        select(symbol, everything()) %>% DataFrame()
    }
  }

  rd = rowData(x)
  rd[,colnames(ann2)] = ann2
  rowData(x) = rd
  return(x)
}


#  **
#  * [description]
#  * @param  {[type]} orgDB [eg: org.Hs.eg.db]
#  * @param  {[character]} gene  the gene name, can be SYMBOL, ENSEMBL, UNIPROT,ALIAS name or mixed of them, note that: when both SYMBOL and ALIAS have a  ENTREZID, will select the SYMBOL mapped ENTREZID]
#  * @return {data.frame table}       [description]
#  *
map_to_entrezid <- function(gene, orgDB = org.Hs.eg.db) {

  try_key = intersect(c("SYMBOL", "ENSEMBL", "UNIPROT", "ALIAS"),columns(orgDB))

  ids_lis = try_key %>% lapply(., function(x){
    suppressMessages(try(mapIds(x = orgDB, keys = gene, keytype = x, column = "ENTREZID"), silent = TRUE))
  })
  names(ids_lis) = paste0("from_", try_key)
  ids_lis = ids_lis[ids_lis %>% sapply(., function(x){class(x) != "try-error"})]
  my_ids <- ids_lis
  my_ids1 <- as.data.frame(do.call(cbind, my_ids), stringsAsFactors = F)

  if(length(ids_lis) < 1){
    message("No gene matched in the this giving organism annotation database! Please check input.")
  }else{
    # my_ids <- ids_lis
    # my_ids1 <- as.data.frame(do.call(cbind, my_ids), stringsAsFactors = F)
    from_cols = my_ids1 %>% apply(., 1, function(x){which(!is.na(x))[1] } ) %>% table(.,useNA =("ifany"))
    names(from_cols) = colnames(my_ids1)[as.numeric(names(from_cols))]
    from_cols1 = from_cols[!(is.na(names(from_cols)))]
    for(i in 1:length(from_cols1)){
      cat(paste0(from_cols1[i]," gene(s) transfrom to ENTREZID ", names(from_cols1)[i] %>% gsub("_"," ",.), ". "))
    }
    cat("\n")
    if(any(is.na(names(from_cols))))
      cat(from_cols[which(is.na(names(from_cols)))] ,"gene(s) do not match in database.\n")
  }

  my_ids1$id <- apply(my_ids1, 1, function(i){
    if(all(is.na(i))) {id = NA
    } else {
      id =  i[which(!is.na(i))[1]]
    }
    return(id)
  })

  return(my_ids1)
}


#' ORA for differenatial test result
#'
#' Enrich biological functions on significant candidate via a over representation analysis.
#'
#' @param x A SummarizedExperiment/DEGdata output from add_adjections or
#' a charachter vector containing candidate identifier(SYMBOL, EntrezID, UniprotID or ENSEMBL).
#' @param contrasts Character, analyse results in which contrasts.
#' @param type Character, one of "GO","KEGG","REACTOME". The datasets for enrichment analysis.
#' @param species The species name.
#' @param by_contrast Logical(1). If true, draw enrichment on each contrast, else draw on the total significant candidates.
#' @param pAdjustMethod Character, one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param ... Other parameters in \code{\link[clusterProfiler]{enricher}()} except the cutoff setting
#'
#' @return
#' A enrichResult or compareClusterResult object according by_contrast.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example
#' data(Silicosis_pg)
#' data <- Silicosis_pg
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' # Make SummarizedExperiment
#' ecols <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, ecols,mode = "delim")
#'
#' # Filter and normalize
#' filt <- filter_se(se, thr = 0, fraction = 0.4, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
#' norm <- normalize_vsn(filt)
#'
#' # Impute missing values using different functions
#' imputed <- impute(norm, fun = "MinProb", q = 0.05)
#'
#' # Test for differentially expressed proteins
#' diff <- test_diff(imputed, type = "control", control  = c("PBS"), fdr.type = "Storey's qvalue")
#' dep <- add_rejections(diff, alpha = 0.01,lfc = 2)
#'
#' # GO enrichment
#' check_organismDB_depends(organism = "Mouse") # check annotation package of Mouse
#' res_ora <- test_ORA(dep, contrasts = "W4_vs_PBS", species = "Mouse",type = "GO")
#' enrichplot::dotplot(res_ora)
#' }
test_ORA <- function(x,
                     contrasts = NULL,
                     type = c(
                       "GO", "KEGG", "REACTOME"),
                     species = "Human",
                     by_contrast = FALSE,
                     pAdjustMethod = c("BH","holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                     ...
){
  check_pak <- check_enrichment_depends()
  if(!isTRUE(check_pak)){
    stop("Packages ",paste0(check_pak,collapse = ", ")," are required for test_ORA, but not found")
  }


  assertthat::assert_that(class(x) == "SummarizedExperiment"|class(x) == "DEGdata",
                          is.null(contrasts)|is.character(contrasts),
                          is.character(species) && length(species) ==1,
                          is.logical(by_contrast) && length(by_contrast) == 1,
                          is.character(pAdjustMethod))
  pAdjustMethod = match.arg(pAdjustMethod)
  type = match.arg(type)


  the_annoSpecies_df = annoSpecies_df()

  if(is.null(contrasts)){
    contrasts = get_contrast(x)
  }
  if(by_contrast && length(contrasts) == 1){
    message("Only contain one contrasts, 'by_contrast' is meanless.")
    by_contrast = F
  }

  if((!by_contrast)){
    if(class(x) == "character"){
      gene_id = x
    }else if(class(x) == "DEGdata" | class(x) == "SummarizedExperiment"){
      # contrasts = get_contrast(x)
      sig = get_signicant(x, return_type = "name")
      # sig = row.names(sig)
    }
    gene = sig ## a character vector
  }else{
    if( class(x) == "character" )
      stop("x should be a test result from add_adjections(), when by_contrast is TRUE")
    # lapply(contrasts, get_signicant, object = x, return_type = "name")
    gene = map(contrasts, get_signicant, object = x, return_type = "name")
    if(length(contrasts) == 1){
      gene = gene[[1]]
      message("only one contrasts exists in", deparse(substitute(x)), ".")
    }
    names(gene) = contrasts ## a list of sig genes in different contrasts
  }

  orgDB = the_annoSpecies_df$pkg[the_annoSpecies_df$species == species]
  require(orgDB, character.only = TRUE)
  orgDB = get(orgDB)
  if(!by_contrast){ ## enricher
    gene_id_table <- map_to_entrezid(gene, orgDB = orgDB)
    ids <- gene_id_table %>% tibble::rownames_to_column() %>%
      dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID) %>%
      filter(!is.na(ENTREZID)) %>% filter(!duplicated(ENTREZID))
    if(type == "GO"){
      cat("Star enrich GO terms by",nrow(ids),"ENTREZIDs.")
      enrich_res <- goAnalysis(ids, organism = species, species_df = the_annoSpecies_df, pAdjustMethod = pAdjustMethod, ...)
    }else if(type == "KEGG"){
      cat("Star enrich KEGG terms by",nrow(ids),"ENTREZIDs.")
      enrich_res <- keggAnalysis(ids, organism = species, species_df = the_annoSpecies_df, pAdjustMethod = pAdjustMethod, ...)
    }else if(type == "REACTOME"){
      cat("Star enrich REACTOME terms by",nrow(ids),"ENTREZIDs.")
      enrich_res <- reactAnalysis(ids, organism = species, species_df = the_annoSpecies_df, pAdjustMethod = pAdjustMethod, ...)
    }

  }else{ ## comparecluster, the input gene_id_list is a list
    gene_id_list <- gene %>% lapply(.,map_to_entrezid, orgDB = orgDB)
    gene_id_list = gene_id_list %>% lapply(., function(gene_id_table){
      gene_id_table %>% tibble::rownames_to_column() %>%
        dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID) %>% filter(!is.na(ENTREZID))
    })
    ## filter out 0 ID vector
    gene_id_list = gene_id_list[which(sapply(gene_id_list,nrow) > 0)]

    if(type == "GO"){
      cat("Star enrich GO terms on",length(gene_id_list),"clusters:", paste0(names(gene_id_list), collapse = ";"),
          "with", paste0(sapply(gene_id_list,nrow), collapse = ";"),"ENTREZIDs.\n")
      enrich_res <- goAnalysis(gene_id_list, organism = species, species_df = the_annoSpecies_df, pAdjustMethod = pAdjustMethod)
    }else if(type == "KEGG"){
      cat("Star enrich KEGG terms on",length(gene_id_list),"clusters:", paste0(names(gene_id_list), collapse = ";"),
          "with", paste0(sapply(gene_id_list,nrow), collapse = ";"),"ENTREZIDs.\n")
      enrich_res <- keggAnalysis(gene_id_list, organism = species, species_df = the_annoSpecies_df, pAdjustMethod = pAdjustMethod)
    }else if(type == "REACTOME"){
      cat("Star enrich REACTOME terms on",length(gene_id_list),"clusters:", paste0(names(gene_id_list), collapse = ";"),
          "with", paste0(sapply(gene_id_list,nrow), collapse = ";"),"ENTREZIDs.\n")
      enrich_res <- reactAnalysis(gene_id_list, organism = species, species_df = the_annoSpecies_df, pAdjustMethod = pAdjustMethod)
    }
  }

  return(enrich_res)
}


#' GSEA data
#'
#' Enrich biological functions on significant candidate via a over representation analysis.
#'
#' @param type Character, one of "GO","KEGG","REACTOME" and "MSigDB". The database for enrichment analysis.
#' @param topn Integer(1), only use topn list with most significant foldchange
#' @param category,subcategory Character. Work when \code{type} is "MSigDB". Use which subset of MSigDB.
#' You can run \code{msigdbr::msigdbr_collections()} to get options.
#' @param ... Other parameters in \code{\link[clusterProfiler]{GSEA}()} except the cutoff setting
#' @inheritParams test_ORA
#' @return
#' A gseaResult object of
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
#' # GSEA
#' check_organismDB_depends(organism = "Mouse") # check annotation package of Mouse
#' res_gsea <- test_GSEA(diff, contrasts = "W4_vs_PBS", species = "Mouse",type = "GO")
#' }
test_GSEA <- function(x,
                      type = c("GO", "KEGG", "REACTOME",
                               "MSigDB"),
                      species = "Human",
                      contrasts = NULL,
                      by_contrast = FALSE,
                      topn= NULL,
                      pAdjustMethod = "BH",
                      category = NULL, subcategory = NULL,
                      # pvalueCutoff = 0.05,  qvalueCutoff = 0.2,
                      ...
){
  check_pak <- check_enrichment_depends()
  if(!isTRUE(check_pak)){
    stop("Packages ",paste0(check_pak,collapse = ", ")," are required for test_GSEA, but not found")
  }

  assertthat::assert_that(class(x) == "SummarizedExperiment"|class(x) == "DEGdata",
                          is.null(contrasts)|(is.character(contrasts)&length(contrasts) ==1),
                          is.character(species) && length(species) ==1,
                          is.logical(by_contrast) && length(by_contrast) == 1,
                          is.character(pAdjustMethod))
  type = match.arg(type)

  the_annoSpecies_df = annoSpecies_df()

  ## extract exist contrasts
  if(class(x) == "SummarizedExperiment"|class(x) == "DEGdata"){
    exist_contrasts = get_contrast(x)
    rd = rowData(x)
  }else if(class(x) == "data.frame"){
    exist_contrasts = colnames(x)[-1]
  }

  ## check contrasts and by_contrast.
  if(is.vector(x)){ ## when x is vector
    message(deparse(substitute(x))," is a vector, 'by_contrast' and contrasts is meanless.")
    by_contrast = F
    contrasts = NULL

  }else{ ## when x is S4 or df
    if(!by_contrast){ ## by_contrast is FALSE
      if(is.null(contrasts) || length(contrasts) < 1){
        contrasts = exist_contrasts[1]
        message("No input contrast and by_contrast is FALSE, use the foldchange in the first contrast: ",contrasts)
      }else if(length(contrasts) > 1 ){
        contrasts = intersect(contrasts, exist_contrasts)[1]
        if(length(contrasts) < 1) stop("Input contrast do not exist, please check")
        message("Only use one contrast when by_contrast is FALSE, but length(contrast)>1,
              use the foldchange in the first contrast: ",contrasts)
      }else if(!contrasts %in% exist_contrasts){
        stop("Input contrast do not exist, please check")
      }

    }else{ ## by_contrast is TRUE, contrasts could be plurality.
      if(is.null(contrasts)|| length(contrasts) < 1){ ## Set default contrasts all exist_contrasts
        contrasts = exist_contrasts
      }else if(length(contrasts) > 1){ ## check input contrasts
        if((!all(contrasts %in% exist_contrasts)) && any(contrasts %in% exist_contrasts)){
          message("Contrasts: ", paste0(contrasts[!contrasts%in%exist_contrasts], collapse = ";"),
                  "do not exist. ")
          contrasts = intersect(contrasts, exist_contrasts)
        }
        if(all(!contrasts %in%exist_contrasts)){
          stop("Input contrasts do not exist in", deparse(substitute(x)))
        }
      }

      if(length(contrasts) == 1){
        message("Only contain useful one contrasts ", contrasts, " 'by_contrast' is meanless.")
        by_contrast = F
      }else{
        message("Use exist contrasts: ", paste0(contrasts, collapse = ";"))
      }
    }
  }
  # x = rowData(my_dep_save)$uH2A_vs_Biotin_diff
  # names(x) = rowData(my_dep_save)$name

  orgDB = the_annoSpecies_df$pkg[the_annoSpecies_df$species == species]
  require(orgDB, character.only = TRUE)
  orgDB = get(orgDB)

  if(is.vector(x)){ ## x is the vector of fc, and names is names/IDs
    genelist = as.numeric(x)
    names(genelist) = names(x)
    genelist = genelist[!is.na(genelist)]
    gene_id_table <- map_to_entrezid(names(x), orgDB = orgDB)
    ids <- gene_id_table %>% tibble::rownames_to_column() %>%
      dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID) %>%
      filter(!is.na(ENTREZID)) %>% filter(!duplicated(ENTREZID))
    if( nrow(ids) < 5 || nrow(ids) < (length(names(x))/2) )
      stop("Mapped gene names is too few!, please check")

    ids1 = ids$ENTREZID[match(names(genelist),ids$name)]
    names(genelist) = ifelse(!is.na(ids1), ids1, names(genelist))
    genelist = genelist[!is.na(ids1)]
    genelist = sort(genelist, decreasing = T) ## order gene list

  }else{ ## x is SE/DEGdata or a dataframe
    if(class(x) == "SummarizedExperiment"|class(x) == "DEGdata"){
      diff_df = rd[,c("name",paste0(contrasts,"_diff"))] %>% as.data.frame
      colnames(diff_df) = gsub("_diff$", "", colnames(diff_df))
    }else if(is.data.frame(x)){
      diff_df = x[,c(1,which(colnames(x) %in% contrasts))]
      diff_df[,-1] =  diff_df[,-1, drop = F] %>% apply(., 2, as.numeric)
    }

    gene_id_table <- map_to_entrezid(diff_df[,1], orgDB = orgDB)
    ids <- gene_id_table %>% tibble::rownames_to_column() %>%
      dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID) %>%
      filter(!is.na(ENTREZID)) %>% filter(!duplicated(ENTREZID))
    if( nrow(ids) < 5 || nrow(ids) < (nrow(diff_df)/2) )
      stop("Mapped gene names is too few!, please check")

    ids1 = ids$ENTREZID[match(diff_df[,1],ids$name)]
    diff_df[,1] = ifelse(!is.na(ids1), ids1, diff_df[,1])
    diff_df = diff_df[!is.na(ids1),,]
    genelist = as.list(diff_df[-1]) %>% lapply(., function(thelist){
      names(thelist) = diff_df[,1]
      thelist = sort(thelist, decreasing = T)
      return(thelist)
    })

    if(length(genelist) == 1)
      genelist = genelist[[1]]
  }

  if(type == "GO"){
    enrich_res <- gsegoAnalysis(genelist, organism = species, pAdjustMethod = pAdjustMethod, ids_table = ids,
                             eps = 1e-12)

  }else if(type == "KEGG"){
    enrich_res <- gsekeggAnalysis(genelist, organism = species, pAdjustMethod = pAdjustMethod, ids_table = ids,
                                  eps = 1e-12)
  }else if(type == "REACTOME"){
    enrich_res <- gsereactAnalysis(genelist, organism = species, pAdjustMethod = pAdjustMethod, ids_table = ids,
                                   eps = 1e-12)
  }else if(type == "MSigDB"){
    enrich_res <- gsemsigdbAnalysis(genelist[[1]], organism = species,
                                    category = category, subcategory = subcategory,
                                    pAdjustMethod = pAdjustMethod, ids_table = ids,
                                    eps = 1e-12)
  }


  return(enrich_res)
}


#' Extract significant enrichment terms base on giving threshold
#'
#' \code{get_ORA_result} filter the enrichment result
#' from \code{test_ORA} though certain threshold.
#'
#' @param reat The output from test_ORA
#' @param ont One of "ALL", "BP", "MF", "CC"
#' @param pvalueCutoff Numeric(1), the p.value cutoff on enrichment result
#' @param qvalueCutoff Numeric(1), the qvalue cutoff on enrichment tests
#' @param simplify Logical(1), if simplify GO terms by \code{\link[clusterProfiler]{simplify}}
#' @param simplify.cutoff Numeric(1), the cutoff value transmitted to  \code{\link[clusterProfiler]{simplify}}
#' @param return_table Logical(1), if true return a enrichResult or a result table
#' @return
#' A enrichResult/compareClusterResult of significant enrichment,
#' or a result table of the significant enrichment if return_table is TRUE.
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
#'
#' res_ora <- test_ORA(dep, contrasts = "W4_vs_PBS", species = "Mouse",type = "GO")
#' res_ora2 <- get_ORA_result(res_ora)
#' }
get_ORA_result <- function(ORA_enrichment, ont = NULL,
                           pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                           simplify = FALSE, simplify.cutoff = 0.7,
                           simplify.measure = c("Wang","Resnik", "Lin", "Rel", "Jiang"),
                           simplify.semData = NULL,
                           return_table = F){
  assertthat::assert_that(class(ORA_enrichment) %in% c("enrichResult", "gseaResult", "compareClusterResult"),
                          is.null(ont) || ont%in% c("GOALL","BP","CC","MF"),
                          is.numeric(pvalueCutoff) && length(pvalueCutoff) == 1,
                          is.numeric(qvalueCutoff) && length(qvalueCutoff) == 1,
                          is.numeric(simplify.cutoff) && length(simplify.cutoff) == 1,
                          is.logical(simplify) && length(pvalueCutoff) == 1,
                          is.logical(return_table) && length(return_table) == 1,
                          is.character(simplify.measure),
                          is.null(simplify.semData) || class(simplify.semData) == "GOSemSimDATA "
                          )
  simplify.measure = match.arg(simplify.measure)
  enrich_type = "other"
  if(class(ORA_enrichment) == "enrichResult"){
    if(ORA_enrichment@ontology %in% c("GOALL","BP","CC","MF")){
      enrich_type ="GO"
    } else if(!is.null(ont)){
      message("'ont' only works on GO enrichResult.")
      ont = NULL
    }
  }else if(class(ORA_enrichment) %in% c("compareClusterResult", "gseaResult")){
    if( all((ORA_enrichment@compareClusterResult$ONTOLOGY %>% unique) %in% c("BP","CC","MF"))){
      enrich_type ="GO"
      exist_onts = ORA_enrichment@compareClusterResult$ONTOLOGY %>% unique()
    } else if(!is.null(ont)){
      message("'ont' only works on GO enrichResult.")
      ont = NULL
    }
  }

  filter_enrichResult <- function(x, pvalueCutoff = 1, qvalueCutoff = 1, ont = NULL){
    if(!(is.null(ont) || ont == "GOALL" || ont == "ALL")){
      x = clusterProfiler.dplyr::filter(x, ONTOLOGY == ont)
      if("ontology" %in% slotNames(ORA_enrichment))
        x@ontology = ont
      if(class(x) == "compareClusterResult")
        x@.call$ont = ont
    }
    x = clusterProfiler.dplyr::filter(x, pvalue <= pvalueCutoff, qvalue <= qvalueCutoff)
    if(pvalueCutoff %in% slotNames(ORA_enrichment))
      x@pvalueCutoff = pvalueCutoff
    if(qvalueCutoff %in% slotNames(ORA_enrichment))
      x@qvalueCutoff = qvalueCutoff
    return(x)
  }
  sig_res <- filter_enrichResult(ORA_enrichment, pvalueCutoff, qvalueCutoff, ont)
  if(nrow(sig_res) == 0){
    message("No result under threashold.")
    return(sig_res)
  }

  if(simplify) {
    if(enrich_type == "GO"){
      if(class(sig_res) == "enrichResult"){
        exist_onts = sig_res@result$ONTOLOGY  %>% unique()
        sig_list <- exist_onts %>%
          sapply(., filter_enrichResult,
                 x= sig_res, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
        sig_list <- lapply(sig_list, clusterProfiler::simplify,
                           cutoff = simplify.cutoff, measure = simplify.measure, semData = simplify.semData)
        sig_list <- lapply(sig_list, function(x){x@result})
        sig_res@result = do.call(rbind,sig_list)

      }else if(class(sig_res) == "compareClusterResult"){
        exist_onts = sig_res@compareClusterResult$ONTOLOGY %>% unique()
        sig_list <- exist_onts %>%
          sapply(., filter_enrichResult,
                 x= sig_res, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
        sig_list <- lapply(sig_list, clusterProfiler::simplify,
                           cutoff = simplify.cutoff, measure = simplify.measure, semData = simplify.semData)
        sig_list <- lapply(sig_list, function(x){x@compareClusterResult})
        sig_res@compareClusterResult = do.call(rbind,sig_list)
      }
    }else{
      message("Simplify function from clusterProfiler package only works on GO enrichResult, skip simplify on other enrichResult")
    }
  }

  if(return_table)
    return(as.data.frame(sig_res))
  return(sig_res)
}



# temp = get_ORA_result(ORA_cluster_res, simplify = T)

# GO_result_cutoff <- function(enrich_result,)

set_readable <- function(ids_vector, ids_table, sep = "/"){
  sapply(ids_vector, function(x){
    IDS = strsplit(x, sep)[[1]]
    ids_table$name[match(IDS,ids_table$ENTREZID)] %>% paste0(., collapse =sep)
  })
}


## ORA functions ----

goAnalysis <- function( gene_id, organism="Human", species_df,
                        pvalueCutoff = 1, qvalueCutoff = 1,
                        ...){
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  if(class(gene_id) == "data.frame"){
    reat_ALL <- try(clusterProfiler::enrichGO(gene = gene_id$ENTREZID, OrgDb = orgDB, ont = "ALL",
                                              pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = F,
                                              ...))
    if(class(reat_ALL) == "try-error") return(reat_ALL)
    reat_ALL@result$geneID %<>% set_readable(., ids_table = gene_id)
  }else if(class(gene_id) == "list"){
    gene_id_list = gene_id
    gene_clusters <- gene_id_list %>% lapply(., function(x){
      set_names(x$ENTREZID, x$name)
    })
    reat_ALL <- try(clusterProfiler::compareCluster(gene_clusters, fun = "enrichGO",
                                                    OrgDb = orgDB, ont = "ALL",
                                                    pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = F))
    if(class(reat_ALL) == "try-error") return(reat_ALL)
    gene_id = do.call(rbind, gene_id_list) %>% .[!duplicated(.$ENTREZID),]
    reat_ALL@compareClusterResult$geneID %<>% set_readable(., ids_table = gene_id)
  }

  reat_ALL@readable = TRUE
  reat_ALL@gene2Symbol = set_names(gene_id$name, gene_id$ENTREZID )
  return(reat_ALL)
}


keggAnalysis <- function(gene_id,
                         organism="Human",
                         pAdjustMethod,species_df,
                         pvalueCutoff = 1, qvalueCutoff = 1,
                         ...){
  organism <- species_df$organism[species_df$species == organism]
  if(class(gene_id) == "data.frame"){
    reat <- try(enrichKEGG(gene = gene_id$ENTREZID, organism = organism,
                           pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                           readable = F,
                           ...), silent = T)
    if(class(reat) == "try-error") return(reat)
    reat@result$geneID %<>% set_readable(., ids_table = gene_id)
  }else if(class(gene_id) == "list"){
    gene_id_list = gene_id
    gene_clusters <- gene_id_list %>% lapply(., function(x){
      set_names(x$ENTREZID, x$name)
    })
    reat <- try(clusterProfiler::compareCluster(geneClusters = gene_clusters, fun = "enrichKEGG",
                               OrgDb = organism, pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                               qvalueCutoff = qvalueCutoff, readable = F,
                               ...), silent = T)
    if(class(reat) == "try-error") return(reat)
    gene_id = do.call(rbind, gene_id_list) %>% .[!duplicated(.$ENTREZID),]
    reat@compareClusterResult$geneID %<>% set_readable(., ids_table = gene_id)
  }

  reat@readable = TRUE
  reat@gene2Symbol = set_names(gene_id$name, gene_id$ENTREZID )
  return(reat)
}



# @import ReactomePA
reactAnalysis <- function(gene_id, organism="Human", pAdjustMethod, species_df,
                          pvalueCutoff = 1, qvalueCutoff = 1, ...){
  organism <- species_df$reactome_organism[species_df$species == organism]
  if(class(gene_id) == "data.frame"){
    reat <- try(ReactomePA::enrichPathway(gene = gene_id$ENTREZID, organism = organism,
                                      pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                                      readable = F,
                                      ...), silent = T)
    if(class(reat) == "try-error") return(reat)
    reat@result$geneID %<>% set_readable(., ids_table = gene_id)

  }else if(class(gene_id) == "list"){
    gene_id_list = gene_id
    gene_clusters <- gene_id_list %>% lapply(., function(x){
      set_names(x$ENTREZID, x$name)
    })
    reat <- try(clusterProfiler::compareCluster(geneClusters = gene_clusters, fun = "enrichPathway",
                                                organism = organism,
                                                pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                                                readable = F,
                                                ...), silent = T)
    if(class(reat) == "try-error") return(reat)
    gene_id = do.call(rbind, gene_id_list) %>% .[!duplicated(.$ENTREZID),]
    reat@compareClusterResult$geneID %<>% set_readable(., ids_table = gene_id)
  }

  reat@readable = TRUE
  reat@gene2Symbol = set_names(gene_id$name, gene_id$ENTREZID )
  return(reat)
}

## GSEA functions ----
# gene_list, a vector. compareCluster worker on gseGO on clusterProfiler(4.3.1.900)
gsegoAnalysis <- function(
  gene_list,
  organism="Human",
  # species_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  ids_table,
  minGSSize = 5, verbose = F, seed = F, eps = 1e-10,
  ...
  ){

  species_df = annoSpecies_df()
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  if(!is.list(gene_list) && is.vector(gene_list)){
    reat_ALL <- try(clusterProfiler::gseGO(gene = gene_list, OrgDb = orgDB, ont = "ALL",
                          pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                          minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                          ...), silent = TRUE)
    if(class(reat_ALL) == "try-error")
      return(reat_ALL)
    reat_ALL@result$core_enrichment %<>% set_readable(., ids_table = ids_table)
  }else if(is.list(gene_list)){
    cat("wdw")
    reat_ALL <- try(clusterProfiler::compareCluster(gene_list, fun = "gseGO",  OrgDb = orgDB, ont = "ALL",
                                                    pAdjustMethod = pAdjustMethod, pvalueCutoff = 1,
                                                    minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                                                    ...
    ), silent = TRUE)
    if(class(reat_ALL) == "try-error")
      return(reat_ALL)
    reat_ALL@compareClusterResult$core_enrichment %<>% set_readable(., ids_table = ids_table)
  }

  reat_ALL@readable = TRUE
  reat_ALL@gene2Symbol = set_names(ids_table$name, ids_table$ENTREZID)
  return(reat_ALL)
}


gsereactAnalysis <- function(gene_list,
                             organism="Human",
                             # species_df,
                             pAdjustMethod = "BH", pvalueCutoff = 1,
                             ids_table,
                             minGSSize = 5, verbose = F, seed = F, eps = 1e-10,
                             ...){
  species_df = annoSpecies_df()
  pkg = species_df$pkg[species_df$species == organism]
  require(pkg, character.only = TRUE)
  orgDB <- get(pkg)

  organism <- species_df$reactome_organism[species_df$species == organism]


  if(!is.list(gene_list) && is.vector(gene_list)){
    reat <- try(ReactomePA::gsePathway(gene = gene_list, organism  = organism,
                           pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                           minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                           ...), silent = TRUE)
    if(class(reat) == "try-error")
      return(reat)
    reat@result$core_enrichment %<>% set_readable(., ids_table = ids_table)
  }else if(is.list(gene_list)){
    reat <- try(clusterProfiler::compareCluster(geneClusters = gene_list, fun = "gsePathway",
                                                organism  = organism,
                                                pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                                                minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                                                ...), silent = TRUE)
    if(class(reat) == "try-error")
      return(reat)
    reat@compareClusterResult$core_enrichment %<>% set_readable(., ids_table = ids_table)
  }

  reat@readable = TRUE
  reat@gene2Symbol = set_names(ids_table$name, ids_table$ENTREZID)
  return(reat)
}


gsekeggAnalysis <- function(gene_list, organism="Human",
                            pAdjustMethod = "BH", pvalueCutoff = 1,
                            # species_df,
                            ids_table,
                            minGSSize = 5, verbose = F, seed = F, eps = 1e-10,
                            ...){
  species_df = annoSpecies_df()
  # pkg = species_df$pkg[species_df$species == organism]
  # require(pkg, character.only = TRUE)
  # orgDB <- get(pkg)
  # cat(pkg)

  organism <- species_df$organism[species_df$species == organism]
  cat(organism)
  if(!is.list(gene_list) && is.vector(gene_list)){
    cat("bbb")
    reat <- try(clusterProfiler::gseKEGG(gene = gene_list, organism = organism,
                                         pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                                         minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                                         ...), silent = TRUE)
    if(class(reat) == "try-error")
      return(reat)
    reat@result$core_enrichment %<>% set_readable(., ids_table = ids_table)
  }else if(is.list(gene_list)){
    cat("aaa")
    reat <- try(clusterProfiler::compareCluster(geneClusters = gene_list, fun = "gseKEGG",
                                                organism  = organism,
                                                pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                                                minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                                                ...), silent = TRUE)
    if(class(reat) == "try-error")
      return(reat)
    reat@compareClusterResult$core_enrichment %<>% set_readable(., ids_table = ids_table)
  }

  reat@readable = TRUE
  reat@gene2Symbol = set_names(ids_table$name, ids_table$ENTREZID)
  return(reat)
}


# @import msigdbr
gsemsigdbAnalysis <- function(gene_list, organism="Human",
                              category = NULL, subcategory = NULL,
                              pAdjustMethod = "BH", pvalueCutoff = 1,
                              ids_table,
                              minGSSize = 5, verbose = F, seed = F, eps = 1e-10,
                              ...){
  assertthat::assert_that(is.character(organism), length(organism) == 1)
  species_df = annoSpecies_df()
  speciesID = species_df$msigdbr_species[species_df$species == organism]
  # cat(speciesID)
  if(!speciesID %in% msigdbr::msigdbr_species()$species_name)
    stop("Organism wrong, please check, it should be one of ",
         paste0(species_df$species[species_df$msigdbr_species != ""], collapse = ", "))

  if(!is.null(category) && !category %in% (msigdbr::msigdbr_collections()$gs_cat %>% unique()))
    stop("category should be NULL or in ",
         paste0((msigdbr::msigdbr_collections()$gs_cat %>% unique()), collapse = ", "))

  if(!is.null(subcategory) && subcategory %in% (msigdbr::msigdbr_collections()$gs_subcat %>% unique()))
    stop("subcategory should be NULL or '' or in ",
         paste((msigdbr::msigdbr_collections()$gs_subcat %>% unique())[-1], collapse = ", "))

  # msigdbr::msigdbr_show_species() %in% species_df$species
  if(!is.null(category)|| (is.null(category) && is.null(subcategory))){
    m_t = msigdbr(species = speciesID, category = category, subcategory = subcategory) %>%
      dplyr::rename(category = gs_cat, subcategory = gs_subcat)
    m_t2g = m_t %>% dplyr::select(gs_name, entrez_gene)
  }else if(!is.null(subcategory)){
    m_t = msigdbr(species = speciesID, category = NULL, subcategory = NULL) %>%
      filter(
        # gs_cat %in% subcategory |
          gs_subcat %in% subcategory) %>%
      dplyr::rename(category = gs_cat, subcategory = gs_subcat)
  }

  m_t2g = m_t %>% dplyr::select(gs_name, entrez_gene)
  if(nrow(m_t2g) == 0) stop("No such dataset of this sepcies or link failure.")


  if(!is.list(gene_list) && is.vector(gene_list)){
    cat("bbb")
    reat <- try(clusterProfiler::GSEA(gene = gene_list, TERM2GENE = m_t2g,
                                      pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                                      minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                                      ...), silent = TRUE)
    if(class(reat) == "try-error")
      return(reat)
    reat@result$core_enrichment %<>% set_readable(., ids_table = ids_table)
    m_t2 = m_t[match(reat@compareClusterResult$Description, m_t$gs_name),] %>% select(category,subcategory)
    reat@compareClusterResult %<>% cbind(m_t2)
  }else if(is.list(gene_list)){
    cat("aaa")
    reat <- try(clusterProfiler::compareCluster(geneClusters = gene_list, fun = "GSEA",
                                                TERM2GENE = m_t2g,
                                                pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff,
                                                minGSSize = minGSSize, verbose = verbose, seed = seed, eps = eps,
                                                ...), silent = TRUE)
    if(class(reat) == "try-error")
      return(reat)
    reat@compareClusterResult$core_enrichment %<>% set_readable(., ids_table = ids_table)
    m_t2 = m_t[match(reat@compareClusterResult$Description, m_t$gs_name),] %>% select(category,subcategory)
    reat@compareClusterResult %<>% cbind(m_t2)
  }

  reat@readable = TRUE
  reat@gene2Symbol = set_names(ids_table$name, ids_table$ENTREZID)
  return(reat)
}





