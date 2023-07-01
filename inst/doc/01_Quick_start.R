## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(BiocStyle)

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(DEP2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

## ----example data-------------------------------------------------------------
## The proteomics data represented in protein-level format (proteingroup.txt output from MaxQuant).
data(Silicosis_pg) 

## The proteomics data represented in peptide-level format (peptides.txt output from MaxQuant).
data(Silicosis_peptide) 

## Experiment design table providing annotations for the samples.
data(Silicosis_ExpDesign)

## Phosphoproteomics data (Phospho(STY)Sites.txt output from MaxQuant).
data(Silicosis_phos) 

## RNA-seq counts data.
data(Silicosis_counts)

## ----start_PG,message=FALSE---------------------------------------------------
## 1. Construct SummarizedExperiement
unique_pg <- make_unique(Silicosis_pg, names = "Gene.names", ids = "Protein.IDs")
ecols <- grep("LFQ.intensity.", colnames(unique_pg))
se <- make_se(unique_pg,                       
              columns = ecols, 
              expdesign = Silicosis_ExpDesign
             )

## 2. Filter
filt <- filter_se(se,
                  thr = 1,
                  filter_formula = ~ Reverse != '+' & Potential.contaminant !="+"
                  )

## 3. Normalize
norm <- normalize_vsn(filt)

## 4. Impute
filt <- DEP2::impute(norm,fun = "MinDet")

## 5. Differential test
diff <- test_diff(filt,type = "control", control = "PBS")

dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

## ----start_PTM,message=FALSE--------------------------------------------------
## 1. Construct SummarizedExperiement
unique_ptm <- make_unique_ptm(Silicosis_phos, 
                              gene_name = "Gene.names", protein_ID = "Protein",
                              aa = "Amino.acid", pos = "Position")
ecols <- grep("Intensity.", colnames(unique_ptm))
se_ptm <- make_se(unique_ptm,                       
                  columns = ecols, 
                  expdesign = Silicosis_ExpDesign
                  )

## 2. Filter
filt_ptm <- filter_se(se_ptm,
                     thr = 1,
                     filter_formula = ~ Reverse!="+" & 
                       Potential.contaminant!="+" & 
                       Localization.prob > 0.7
                     )

## 3. Normalize
norm_ptm <- normalize_vsn(filt_ptm)

## 4. Impute
filt_ptm <- DEP2::impute(norm_ptm, fun = "QRILC")

## 5. Differential test
diff_ptm <- test_diff(filt_ptm,type = "control", control = "PBS")

dep_ptm <- add_rejections(diff_ptm, alpha = 0.01)

## ----start_pep, eval=F--------------------------------------------------------
#  ## 1. Construct SummarizedExperiement
#  ecols <- grep("Intensity.", colnames(Silicosis_peptide), value = TRUE)
#  pe = make_pe_parse(Silicosis_peptide,
#                      columns = ecols,   # columns is the 'Intensity' colunmns
#                      mode = "delim", sep = "_",
#                      remove_prefix = TRUE
#                     )
#  
#  ## 2. Filter
#  pe = filter_pe(pe,
#                 thr = 1,
#                 fraction = 0.3,
#                 filter_formula = ~ Reverse != '+' & Potential.contaminant !="+"
#                 )
#  
#  ## 3. Impute
#  pe <- impute_pe(pe, fun = "QRILC", name = "peptideImp")
#  
#  ## 4. Normalize
#  pe <- normalize_pe(pe,method = "quantiles.robust", i = "peptideImp")
#  
#  ## 5. Aggregate peptides quantity to proteins'
#  begin_time = Sys.time()
#  pe <- aggregate_pe(pe, fcol = "Proteins", reserve = "Gene.names")
#  print(Sys.time() - begin_time) ## Required few minutes
#  
#  ## 6. Transform a SummarizedExperiement of protein quantities.
#  se <- pe2se(pe, names = "Gene.names", ids = "smallestProteingroups")
#  
#  
#  ## 7. Differential test
#  diff_pep <- test_diff(se,type = "control", control = "PBS")
#  dep_pep <- add_rejections(diff_pep, alpha = 0.01)

## ----start_rna,message=FALSE--------------------------------------------------
## 0. Check the depend packages
check_RNAseq_depends()

## 1. Construct DESeqDataSet 
dds <- make_dds_parse(Silicosis_counts,mode = "delim")

## 2. Filter
dds <- filter_se(dds, fraction = 0.3, thr = 1, rowsum_threshold = 35)

## 3. Transform ID (optional)
dds <- ID_transform(dds, fromtype = "ENSEMBL", species = "Mouse")
head(rownames(dds))

## 3. Differential test
diff <- test_diff_deg(dds, type = "control", control = "PBS")
deg <- add_rejections(diff, alpha = 0.01, lfc = 1)

## ----heatmap1-----------------------------------------------------------------
class(dep)

## Sort the conditions
dep = Order_cols(dep, order = c("PBS","W2","W4","W6","W9","W10"))

## Heatmap on centered values
DEP2::plot_heatmap(dep,
                   cluster_column = F,
                   kmeans = T,        # kmeans cluster rows
                   col_limit = 4      # limit of the color scale
                   )     

## Heatmap on contrast log2 foldchange
DEP2::plot_heatmap(dep,
                   type = "contrast",
                   cluster_column = F,
                   kmeans = T, 
                   col_limit = 4)

## Manual contrast
DEP2::plot_heatmap(dep,
                   manual_contrast = c("W6_vs_PBS","W4_vs_PBS"),
                   cluster_column = F,
                   kmeans = T,        
                   col_limit = 4)

## Change color panel 
DEP2::plot_heatmap(dep,
                   type = "contrast",
                   color = "PiYG",
                   cluster_column = F,
                   kmeans = T,
                   col_limit = 4)


## Parameters in ComplexHeatmap::Heatmap are available
ht <- DEP2::plot_heatmap(dep,
                         cluster_column = F,
                         kmeans = T,
                         cluster_rows = F
                         )

## Heatmap can be saved as follow
# pdf("ht.pdf")
# ht
# dev.off()

## ----volcano2-----------------------------------------------------------------
plot_volcano(dep)


plot_volcano(deg, 
             contrast = "W6_vs_PBS",  # contrast
             label_number = 20
             )

library(ggplot2)
(volc <- plot_volcano(deg, contrast = "W4_vs_PBS") + xlim(-10,10))

# ggplot2::ggsave(filename = "volcano.png") # export to file

## ----export_table-------------------------------------------------------------
## Significant table
sig_tb <- get_signicant(dep,
                        return_type = "table" 
                        )
sig_tb_W6_vs_PBS <- get_signicant(dep,
                                  return_type = "table" ,
                                  contrasts = "W6_vs_PBS"
                                  )

## Full table
full_tb <- get_df_wide(dep)

dim(sig_tb)
dim(sig_tb_W6_vs_PBS)
dim(full_tb)

## Save table in R
# write.csv(full_tb,"fulldataset.cs",row.names = F)

## ----sessioninfo, echo=FALSE--------------------------------------------------
sessionInfo()

