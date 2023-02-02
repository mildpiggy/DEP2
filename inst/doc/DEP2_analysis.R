## ----opts ,include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = T,
  eval = T,
  comment = "#>"
  # ,fig.path = "man/figures/README-"
  # ,out.width = "100%"
)

## ----library, include=FALSE---------------------------------------------------
suppressPackageStartupMessages({
  # devtools::load_all()
  library(BiocStyle)
  library(tidyr)
  library(DEP2)
  library(SummarizedExperiment)
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(missForest)
  library(patchwork)
})
# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

## ----readin,collapse=T--------------------------------------------------------
## a proteinGroups data of a silicosis mouse model
data(Silicosis_pg) 
colnames(Silicosis_pg)

## Formatting name(gene symbol) and id(protein ID). 
## Generate a unique names for each protein. 
unique_pg <- make_unique(Silicosis_pg, names = "Gene.names", 
                         ids = "Protein.IDs", delim = ";") # names and ids are columns in table
## The duplicated names are uniqued by increasing decimal number.
unique_pg %>% filter(Gene.names!="") %>% group_by(Gene.names) %>% summarize(frequency = n(), names = paste0(name,collapse = ";")) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1) %>% head

## ----SE,collapse=T------------------------------------------------------------
## Take expression columns(LFQ intensity in this cases).
ecols <- grep("LFQ.intensity.", colnames(unique_pg))

## Construct SE
### 1. import a experiment design table.
data(Silicosis_ExpDesign)
expdesign <- Silicosis_ExpDesign
# expdesign <- read.delim(system.file("extdata/Silicosis_expdesign.txt", package="DEP2"))
head(expdesign,6) # the example experiment design

se_pg <- DEP2::make_se(unique_pg, columns = ecols, expdesign = expdesign, 
                       log2transform = T)
class(se_pg)
head(assay(se_pg), 5) ## log2-transformed expression assay 

### 2. Generate a experiment design by parsing expression column names.
se_pg2 <- DEP2::make_se_parse(unique_pg, columns = ecols, mode = "delim", 
                              sep = "_", remove_prefix = T, log2transform = T)

## ----readinpep,collapse=T-----------------------------------------------------
## a peptides table from a silicosis mouse model.
data(Silicosis_peptide) 

## Here we just select two groups to save computational cost in these vignette
(ecols <- grep("Intensity.PBS|Intensity.W4", colnames(Silicosis_peptide), value = T)) 
pe_peptides <- make_pe_parse(Silicosis_peptide, columns = ecols,   # columns is the abundance columns
                             mode = "delim", sep = "_",
                             remove_prefix = T, log2transform = T) # log2transform 
pe_peptides # a QFeatures object, with a peptideRaw assay
colData(pe_peptides)

## ----readinpho----------------------------------------------------------------
## phosphorylated peptides table of the silicosis mouse model.
data(Silicosis_phos) 

## Format the modification information and generated modified-peptides identifier.
## aa and pos is the modified amino acids and modified site in protein.
unique_pho <- make_unique_ptm(Silicosis_phos, gene_name = "Gene.names", 
                              protein_ID = "Protein", aa = "Amino.acid",
                              pos = "Position") 

## ----readinpho2---------------------------------------------------------------
## The formatted name and ID for PTM
unique_pho %>% select("name","ID","gene_name","protein_ID","modified_aa","modified_pos") %>% head(7)

## Take expression columns.
ecols <- grep("Intensity.", colnames(unique_pho))

## Construct a SE object mentioned before.
se_ptm <- DEP2::make_se_parse(unique_pho, columns = ecols, mode = "delim", sep = "_", remove_prefix = T, log2transform = T)

## ----filter pg, collapse=T----------------------------------------------------
filter_pg <- filter_se(se_pg,
                       thr = 1,  ## the threshold of missing number in at least one condition
                       fraction = 0.3, ## the threshold of missing occupancy in each protein
                       filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" ## filter upon Reverse and contaminant
                      )
## Further filter can perform though filter_formula.
## For example, keep proteingroups with at least one unique peptide in further analysis.
filter_pg <- filter_se(filter_pg,
                       filter_formula = ~ Unique.peptides > 0)

get_df_wide(se_pg)$Reverse %>% table
get_df_wide(filter_pg)$Reverse %>% table

(plot_coverage(se_pg) + ggtitle("Before filter")) + 
  (plot_coverage(filter_pg) + ggtitle("After filter"))

(plot_frequency(se_pg) + ggtitle("Identification overlap before filter")) / 
  (plot_frequency(filter_pg) + ggtitle("Identification overlap after before filter"))

## ----message=FALSE------------------------------------------------------------
norm_pg <- normalize_vsn(filter_pg)
plot_normalization(filter_pg, norm_pg)

## ----missing heatmap, message=FALSE-------------------------------------------
plot_missval(filter_pg)

## ----impute, collapse=T-------------------------------------------------------
set.seed(35)
sample_rows <- sample(1:nrow(norm_pg), 300)
norm_pg_sample = norm_pg[sample_rows,] # random sample 150 features to reduce runing time
plot_detect(norm_pg_sample)

## Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
imp_pg_MinProb <- DEP2::impute(norm_pg_sample, fun = "MinProb", q = 0.01)

## Impute missing data using k-nearest neighbour approach
imp_pg_knn <- DEP2::impute(norm_pg_sample, fun = "knn")

## Impute missing data using missForest
imp_pg_RF <- DEP2::impute(norm_pg_sample, fun = "RF", ntree = 50, mtry = 5) 

## Impute missing data using Gibbs
imp_pg_GSimp <- DEP2::impute(norm_pg_sample, fun = "GSimp", hi_q = 0.1,
                       iters_each=40, iters_all=8)

## Impute missing data using QRILC
imp_pg_QRILC <- DEP2::impute(norm_pg_sample, fun = "QRILC") 

## ----message=FALSE------------------------------------------------------------
plot_detect(norm_pg_sample)

NAs <- is.na(assay(norm_pg_sample)) 
## the imputed values by different methods.
imps <- list("GSimp" = imp_pg_GSimp, "QRILC" = imp_pg_QRILC, "MinProb" = imp_pg_MinProb, "RF" = imp_pg_RF, "knn" = imp_pg_knn) %>% 
  lapply(function(se){
    x = assay(se) %>% data.frame %>% gather("label", "value") %>% 
      left_join(colData(se)[c("label","condition")],copy = T) %>%
      magrittr::extract(as.vector(NAs),)
  }) %>% data.table::rbindlist(idcol = "method")

## the original normalized values without imputation
nonimps <- assay(norm_pg_sample) %>% data.frame %>% gather("label", "value") %>%
  left_join(colData(norm_pg_sample)[c("label","condition")],copy = T) %>% 
  magrittr::extract(!as.vector(NAs),) %>% mutate(method = "non_impute") %>%
  select(method,everything())

library(ggridges)
ggplot(rbind(imps, nonimps),aes(x = value,y = factor(method,level = unique(method)))) + 
  geom_density_ridges(fill = "#027ad450", scale = 1.2,
                      jittered_points = TRUE,position = position_points_jitter(height = 0),
                      point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.7) +
  ylab("Impute method")+ ylab("Log2 value") + xlim(c(9,39))+
  theme_DEP1()

## ----rmtemp-------------------------------------------------------------------
rm(list = c("norm_pg_sample", "imp_pg_knn","imp_pg_MinProb","imp_pg_RF","imp_pg_GSimp"))
## Impute missing data using q-th quantile for following analyze
imp_pg <- impute(norm_pg, fun = "MinDet" ,q = 0.01)

## ----testdiff-----------------------------------------------------------------
## Test every sample versus PBS control
diff_pg <- test_diff(imp_pg, type = "control", control = "PBS", fdr.type = "BH")
## Test on manul contrasts
diff_pg2 <- test_diff(imp_pg, type = "manual", test  = c("W4_vs_PBS"), fdr.type = "Storey's qvalue")

## ----rejections1--------------------------------------------------------------
## Add significant rejection for features, based on 
dep_pg <- add_rejections(diff_pg, alpha = 0.01, lfc = 2)
## get significant subset
dep_pg_sig <- get_signicant(dep_pg)
nrow(dep_pg_sig)

## ----rejections plot 1--------------------------------------------------------
### volcano plot on contrast "W4_vs_PBS"
plot_volcano(dep_pg, contrast = "W4_vs_PBS", adjusted = F)
## plot the cutoff line
plot_volcano(dep_pg, contrast = "W4_vs_PBS", adjusted = F,
             add_threshold_line = "intersect", pCutoff  = 0.05, fcCutoff = 1)

## ----rejections2,fig.height=7-------------------------------------------------
## thresholdmethod = "curve"
dep_pg_curve <- add_rejections(diff_pg, thresholdmethod = "curve", curvature  = 2, x0_fold = 2)

## the cutoff line on volcano
plot_volcano(dep_pg_curve, contrast = "W4_vs_PBS", add_threshold_line = "curve", curvature = 2, x0_fold = 2) /
  plot_volcano(dep_pg_curve, contrast = "W4_vs_PBS", add_threshold_line = "curve", curvature = 2, x0_fold = 2)


## -----------------------------------------------------------------------------
## Check the fitted Gaussian curve
plot_diff_hist(dep_pg_curve, contrasts = c("W4_vs_PBS", "W6_vs_PBS"))

## Get the fit result 
plot_diff_hist(dep_pg_curve, plot = F) ## a table of gaussian args \sigma and \mu

## ----heatmap,collapse=T-------------------------------------------------------
plot_heatmap(dep_pg)

## Reorder columns by condition
dep_pg = DEP2::Order_cols(dep_pg,order = c("PBS","W2","W4","W6","W9","W10"))
plot_heatmap(dep_pg, cluster_columns = F, kmeans = T, k = 5, seed = 1) # cluster features

## Only plot the clusters that are up-regulated in treatment groups.
plot_heatmap(dep_pg, cluster_columns = F, kmeans = T, k = 5, seed = 1, col_limit = 4,
             split_order = c(1,2,5)
             )
## plot on select contrast
plot_heatmap(dep_pg, manual_contrast = "W4_vs_PBS")

## ----filter pe----------------------------------------------------------------
filt_pe <- filter_pe(pe_peptides, 
                     thr = 1,fraction = 0.4,  ## filter based on missing values
                     filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" )

## ----imp pe, message=FALSE----------------------------------------------------
## Use QFeatures::addAssay() and DEP2::impute
pe <- QFeatures::addAssay(filt_pe,
                              DEP2::impute(filt_pe[["peptideRaw"]], fun = "QRILC"), 
                              name = "peptideImp")
plot_imputation(pe[["peptideRaw"]], pe[["peptideImp"]])

## ----imp before norm2---------------------------------------------------------
pe <- DEP2:::normalize_pe(pe,method = "quantiles.robust",
                               i = "peptideImp",    ## nomalize on which assay
                               name = "peptideNorm" ## output assay name
                               )

plot_normalization(pe[["peptideRaw"]],pe[["peptideImp"]],pe[["peptideNorm"]])


## ----aggrega------------------------------------------------------------------
begin_time = Sys.time()
pe <- DEP2::aggregate_pe(pe,
                         aggrefun = "RobustSummary",
                         aggregate_Peptide_Type = "Unique + Razor",
                         fcol = "Proteins",                  # the protein ID column in input table
                         peptide_assay_name = "peptideNorm", # the assay to aggregate
                         reserve = "Gene.names"              # reserve the Gene.names in row data.
)

print(Sys.time() - begin_time)

pe # protein result is in the protein assay
se_pep <- pe2se(pe, names = "Gene.names", ids = "smallestProteingroups")

## -----------------------------------------------------------------------------
## Differential test
diff_pep <- DEP2::test_diff(se_pep,type = "control", control = "PBS", fdr.type = "BH")
## Add rejections
dep_pep <- add_rejections(diff_pep,alpha = 0.01, lfc = 2)
get_signicant(dep_pep) %>% nrow

## ----DEP visualization--------------------------------------------------------
## Volcano
plot_volcano(dep_pep, contrast = "W4_vs_PBS", add_threshold_line = "intersect") 
## Heatmap
plot_heatmap(dep_pep, kmeans = T,k = 5,col_limit = 6) 

## ----message=FALSE, warning=FALSE---------------------------------------------
## Filter base on both missing occupancy and the localization probability for this site.
filt_ptm <- filter_se(se_ptm, 
                      thr = 1, fraction = 0.3, 
                      filter_formula = ~ Reverse!="+" & Potential.contaminant!="+"&Localization.prob>0.7)

norm_ptm <- normalize_vsn(filt_ptm)

imp_ptm <- impute(filt_ptm, fun= "knn")

diff_ptm <- test_diff(imp_ptm, type = "manual", test = "PBS_vs_W6" , fdr.type = "BH")

dep_ptm <- DEP2::add_rejections(diff_ptm, alpha = 0.05, lfc = 1)
plot_volcano(dep_ptm,adjusted = T, add_threshold_line = "intersect")

## -----------------------------------------------------------------------------
## Get result table
DT::datatable(head(get_results(dep_pg)), options = list(scrollX = T))
## Get full data set with row annotation.
DT::datatable(head(get_results(dep_pg)), options = list(scrollX = T))
## Get full data set of significant features
get_signicant(dep_pg,return_type = "table") %>% select(ends_with("significant")) %>% summary

## -----------------------------------------------------------------------------
# Test on the contrast W4_vs_PBS.
diff_pg <- test_diff(imp_pg, type = "manual", test  = c("W4_vs_PBS"), fdr.type = "Storey's qvalue")
dep_pg <- diff_pg %>% add_rejections(alpha = 0.05,lfc = 2)
diff_pep <- DEP2::test_diff(se_pep,type = "control", control = "PBS", fdr.type = "Storey's qvalue")
dep_pep <- add_rejections(diff_pep,alpha = 0.05, lfc =2)

# Filter out pg quanfied by only one peptides.
dep_pep = filter_se(dep_pep,filter_formula = ~.n  > 2)
dep_pg = filter_se(dep_pg,filter_formula = ~Razor...unique.peptides  > 2)

# Plot venn
plot_multi_venn(omics_list = list(PG= dep_pg, pep=dep_pep))

# If set the proteins remained in both peptides as background
overlap_protein = intersect(rownames(dep_pep), rownames(dep_pg))
plot_multi_venn(omics_list = list(PG= dep_pg, pep=dep_pep), background = overlap_protein)

## ----mulheatmap---------------------------------------------------------------
# Take some proteins may response to interleukin-1 as an example.
IL1_relative_genes <-c("Irg1", "Il1rn", "Saa3", "Zbp1", "Ccl6", 
                       "Serpine1", "Ccl21a", "Pycard", "Irak2", "Vrk2", 
                       "Fn1", "Il1r1", "Irf1", "Ccl9", "Mapk11", "Tank", "Mapk13") 
# Order the samples by time.
dep_pg <- Order_cols(dep_pg, order = c("PBS","W2","W4","W6","W9","W10"))
plot_multi_heatmap(omics_list = list(PG= dep_pg, pep=dep_pep),
                   choose_name = IL1_relative_genes,
                   width = 4)

## ----enrichment, message=FALSE, warning=FALSE---------------------------------
diff_pg <- test_diff(imp_pg, type = "control", control  = c("PBS"), fdr.type = "Storey's qvalue")
dep_pg <- diff_pg %>% add_rejections(alpha = 0.01,lfc = 2)

check_enrichment_depends() # check required packages
require(org.Mm.eg.db) # organism annotation for mouse
require(clusterProfiler) # organism annotation for mouse

## Over representation analysis
# On one contrast
res_ora <- test_ORA(dep_pg, contrasts = "W4_vs_PBS", species = "Mouse",type = "GO")
class(res_ora)

# On each contrasts
res_ora2 <- test_ORA(dep_pg, by_contrast = T, species = "Mouse",type = "GO")
class(res_ora2)

## Gene set enrichment analysis
res_gsea <- DEP2::test_GSEA(dep_pg, by_contrast = F, species = "Mouse",type = "GO")
class(res_gsea)

## ----dotplot------------------------------------------------------------------
## Dotplot function from DOSE package
dotplot(res_ora)

## ----gesaplot-----------------------------------------------------------------
## GSEA plot on one category
gseaplot(res_gsea,geneSetID = "GO:0071347")

## ----cnet---------------------------------------------------------------------
## Cnetplot from enrichplot
res_gsea2 <- res_gsea %>% 
  filter(ID %in% c("GO:0071347","GO:0035329","GO:0032660","GO:0090263","GO:0035456","GO:0038061"))
cnetplot(res_gsea2,foldChange = res_gsea@geneList,node_label = "category",showCategory = 6)

## ----tcheatmap,fig.height=14--------------------------------------------------
## Perform cluster on signicant protein in all contrast
tc_pg <- get_tc_cluster(get_signicant(dep_pg),
               group_order = c("PBS","W2","W4","W6","W9","W10") # reset the group order
              )
## Heatmap
tc_pg$ht

## ----tcheatmap2---------------------------------------------------------------
## On phospholated sites
tc_ptm <- get_tc_cluster(get_signicant(dep_ptm),
                         group_order = c("PBS","W2","W4","W6","W9","W10"), # reset the group order
                         row_font_size = 0,
                         k = 8,
                         heatmap_height = 2
)
tc_ptm$ht

## ----tcres--------------------------------------------------------------------
## The cluster result table
tc_ptm$res %>% head

## ----PPI1,warning=FALSE,eval=TRUE,message=FALSE-------------------------------
# Check required packages for network construction.
check_PPI_depends()
# Load STRING data. If local STRING is missing, PPI_res will download to system.file("PPIdata", "Mouse",package = "DEP2")
load_PPIdata(speciesname = "Mouse") 

# Load STRING data and perform PPI analysis. If local STRING is missing, PPI_res will download to system.file("PPIdata",package = "DEP2")
PPI_res <- test_PPI(
  dep_pg,
  contrasts = "W4_vs_PBS",
  species = "Mouse",
  STRING.version = "11.5",
  score_cutoff = 400
)

# igraph network
PPI_ig <- PPInetwork(PPI_res,returntype = "igraph")
igraph::plot.igraph(PPI_ig)


## ----PPI2,eval=FALSE,fig.height=8---------------------------------------------
#  ## visNetwork network
#  PPInetwork(PPI_res,returntype = "visNetwork")

## ----sessioninfo, echo=FALSE--------------------------------------------------
sessionInfo()

