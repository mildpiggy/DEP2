## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  comment = "#>"
)
library(BiocStyle)

## ----workplace,eval=FALSE-----------------------------------------------------
#  # Download OmicsExample firstly
#  # Change the path to Omics Example
#  example_path = "the/path/to/OmicsExample"
#  knitr::opts_knit$set(
#    root.dir = example_path
#    )
#  
#  library(DEP2)
#  library(dplyr)

## ----workplace0,include=FALSE-------------------------------------------------
example_path = "../example_data/OmicsExample/"
knitr::opts_knit$set(
  root.dir = example_path
  )
library(DEP2)
library(dplyr)

## ----mq-----------------------------------------------------------------------
# Load proteinGroups table
mq_pgfile <- "./A_spikeIn/MaxQuant_output/proteinGroups.txt.gz"

pg <- read.csv(mq_pgfile,sep = "\t")

## ---- mq_uni------------------------------------------------------------------
## Generate unique names and ids
unique_pg <- make_unique(pg, 
                         name = "Gene.names", #gene 
                         ids = "Protein.IDs"   #protein
                         )

## New columns "name" and "id", which is necessary to function make_se
head(unique_pg[,c("name","ID","Gene.names","Protein.IDs")])

## ----mq_se--------------------------------------------------------------------
## Expression columns in the input data
ecols <- grep("LFQ.intensity.", colnames(unique_pg))

## Experiment design table
expDesign_file <- "./A_spikeIn/MaxQuant_output/expdesign.txt.gz"
expdesign <- read.table(expDesign_file, sep = "\t", header = T)
expdesign

## Convert SE with expdesign table
se <- make_se(unique_pg, columns = ecols, expdesign = expdesign)

## Convert SE with expdesign table using parse function
se <- make_se_parse(unique_pg, columns = ecols, mode = "delim")

## ----get_exdesign-------------------------------------------------------------
## How experiment design is extracted from colnames. Split by delim or by character number 
sample_names = c("Quantity.A_1","Quantity.A_2","Quantity.B_1",
                 "Quantity.B_3","Quantity.B_2","Quantity.A_3")
get_exdesign_parse(sample_names,mode = "delim", sep = "_", remove_prefix = TRUE)

sample_names = c("Quantity.A1","Quantity.A2","Quantity.B1",
                 "Quantity.B3","Quantity.B2","Quantity.A3")
get_exdesign_parse(sample_names,mode = "char", chars = 1, remove_prefix = TRUE)

## ----FragP--------------------------------------------------------------------
# Load combined_protein.tsv table
FragP_pgfile <- "./A_spikeIn/FragPipe_output/combined_protein.tsv.gz"

FragP_pg <- read.csv(FragP_pgfile,sep = "\t")

## ----FragP_unique2------------------------------------------------------------
## Generate unique names and ids
unique_FragP_pg <- make_unique(FragP_pg, 
                               name = "Gene",  #gene 
                               ids = "Protein.ID"   #protein
                               )

## new columns name and id, which is necessary to make_se
head(unique_FragP_pg[,c("name", "ID", "Gene", "Protein.ID")])

## ----FragP_se_frg-------------------------------------------------------------
## Expression col in the combined_protein.tsv
ecols <- grep(".MaxLFQ.Intensity", colnames(unique_FragP_pg),value = T)
ecols

## Convert SE with expdesign table
se <- make_se_parse(unique_FragP_pg, 
                    columns = ecols, 
                    mode = "delim",
                    remove_suffix = TRUE    ## remove the col suffix ".MaxLFQ.Intensity"
                    )

## The suffix is removed in SE
colData(se)
colnames(se)

## ----diann_load1--------------------------------------------------------------
# Load combined_protein.tsv table
Diann_pgfile <- "./A_spikeIn/Diann_output/report.pg_matrix.tsv.gz"

Diann_pg <- read.csv(Diann_pgfile,sep = "\t", fileEncoding="latin1")

## ----unique2------------------------------------------------------------------
## Generate unique names and ids
unique_diann_pg <- make_unique(Diann_pg, 
                               name = "Genes",  #gene 
                               ids = "Protein.Group"   #protein
                               )

## New columns ”name“ and "ID", which is necessary to make_se
head(unique_diann_pg[,c("name", "ID", "Genes", "Protein.Group")])

## ----make_se_diann1-----------------------------------------------------------
## Expression col in the DIA-NN report.pg_matrix
ecols <- grep(".raw$", colnames(unique_diann_pg),value = T)
ecols

## Convert SE with expdesign table
se <- make_se_parse(unique_diann_pg, 
                    columns = ecols, 
                    mode = "delim",
                    remove_prefix = TRUE,   ## remove the file prefix 
                    remove_suffix = TRUE    ## remove the col suffix ".raw"
                    )

## The file prefix and suffix is removed in SE
colData(se)
colnames(se)

## ----diann_load2--------------------------------------------------------------
## Load combined_protein.tsv table
Diann_repfile <- "./A_spikeIn/Diann_output/report.tsv.gz"
# Diann_repfile <- "./example_data/OmicsExample/A_spikeIn/Diann_output/report.tsv.gz"
Diann_rep <- read.csv(Diann_repfile,sep = "\t") #
dim(Diann_rep) # 390479 rows, each row is a precursor.

## Filter out proteingroups exceeding the threshold value
Diann_rep = dplyr::filter(Diann_rep, PG.Q.Value < 0.01)


## Table report.tsv stores both Precursor- and PG-level quality
DT::datatable(head(Diann_rep) %>% mutate_if(is.character, utf8::utf8_encode), options = list(scrollX = T,pageLength = 6))

## ----diann_resha--------------------------------------------------------------

Diann_rep_wided = reshape_long2wide(Diann_rep,
                                    sample_col = "File.Name",      # the column labeling sample names
                                    feature_col = "Protein.Group", # PG IDs
                                    expression_col = "PG.MaxLFQ",  # PG quantity. Normalized one is also ok。
                                    shrink_ident_cols  = "Genes",  # Gene names
                                    extend_ident_cols =  "Protein.Q.Value",  # optional, some identification info.
                                    remove_sample_prefix = FALSE,  # remove prefix in sample_col
                                    remove_sample_suffix = FALSE   # remove suffix in sample_col
                                    )


## ----DT_wide------------------------------------------------------------------
DT::datatable(head(Diann_rep_wided,3), options = list(scrollX = T,pageLength = 3))

## ----diann_unique1------------------------------------------------------------
## Generate unique names and ids
unique_diann_pg2 <- make_unique(Diann_rep_wided, 
                                name = "Genes",  #gene 
                                ids = "Protein.Group"   #protein
                                )

## new columns name and id, which is necessary to make_se
head(unique_diann_pg2[,c("name", "ID", "Genes", "Protein.Group")],4)

## ----make_se_diann3-----------------------------------------------------------
## Expression col in the DIA-NN report.pg_matrix
ecols <- grep(".raw$", colnames(unique_diann_pg2),value = T)
ecols # contains Protein.Q.Value columns

ecols = ecols[1:6]

## Convert SE with expdesign table
se <- make_se_parse(unique_diann_pg2, 
                    columns = ecols, 
                    mode = "delim",
                    remove_prefix = TRUE,   ## remove the file prefix 
                    remove_suffix = TRUE    ## remove the col suffix ".raw"
                    )

## The file prefix and suffix is removed in SE
colData(se)
colnames(se)

## ----spe_load1----------------------------------------------------------------
Spe_repfile <- "./A_spikeIn/Spectronaut_output/DIA_MSStats_Report.xls.gz"
Spe_rep <- read.csv(Spe_repfile,sep = "\t") 
dim(Spe_rep) # 318132 rows, each row is a precursor.

## Store both Precursor-, Peptides- and PG-level quality
DT::datatable(head(Spe_rep) %>% mutate_if(is.character, utf8::utf8_encode), options = list(scrollX = T,pageLength = 6))

## ----spe_resha----------------------------------------------------------------
Spe_rep_wided = reshape_long2wide(Spe_rep,
                                  sample_col = "R.FileName",      # the column labeling sample names
                                  feature_col = "PG.ProteinGroups", # PG IDs
                                  expression_col = "PG.Quantity",   # PG quantity.
                                  shrink_ident_cols  = c("PG.Genes","PG.Qvalue"),  # Gene names
                                  remove_sample_prefix = TRUE,    # remove prefix in sample_col
                                  remove_sample_suffix = FALSE  
                                  )


## ----spe_DT_wide--------------------------------------------------------------
DT::datatable(head(Spe_rep_wided), options = list(scrollX = T,pageLength = 6))

## ----spe_unique2--------------------------------------------------------------
## Generate unique names and ids
unique_spe_pg <- make_unique(Spe_rep_wided, 
                             name = "PG.Genes",  #gene 
                             ids = "PG.ProteinGroups"   #protein
                             )

## new columns name and id, which is necessary to make_se
head(unique_spe_pg[,c("name", "ID", "PG.Genes", "PG.ProteinGroups")],4)

## ----make_se_diann2-----------------------------------------------------------
## Expression columns in the Spectronaut output
ecols <- 2:7
colnames(unique_spe_pg)[ecols]

## Convert SE with expdesign table
se <- make_se_parse(unique_spe_pg, 
                    columns = ecols, 
                    mode = "delim",
                    remove_prefix = TRUE,   ## remove the file prefix 
                    remove_suffix = FALSE   ## no suffix in this case
                    )

## The file prefix and suffix is removed in SE
colData(se)
colnames(se)

## ----MQ_pep-------------------------------------------------------------------
mq_pepfile <- "./A_spikeIn/MaxQuant_output/peptides.txt.gz"

mq_pep <- read.csv(mq_pepfile,sep = "\t")

## ----MQ_pep2------------------------------------------------------------------
ecols <- grep("Intensity.", colnames(mq_pep)) ## the peptides intensity cols
mq_pe <- make_pe(mq_pep, columns = ecols, expdesign = expdesign)
mq_pe

## ----Fragp_pep1---------------------------------------------------------------
FragP_pepfile <- "./A_spikeIn/FragPipe_output/combined_peptide.tsv.gz"

FragP_pep <- read.csv(FragP_pepfile,sep = "\t")

## ----Frag_pe1-----------------------------------------------------------------
ecols <- grep("[0-9].Intensity", colnames(FragP_pep),value = T)  ## the peptides intensity cols

FragP_pe <- make_pe_parse(FragP_pep, columns = ecols, 
                          mode = "delim", 
                          remove_suffix = T ## remove suffix
                          )
colData(FragP_pe)

## ----diann_pep1---------------------------------------------------------------

## filter out proteingroups exceed the threshold value
Diann_rep = dplyr::filter(Diann_rep, PG.Q.Value < 0.01)

Diann_pep_wided = reshape_long2wide(Diann_rep,
                                    sample_col = "File.Name",      # the column labeling sample names
                                    feature_col = "Stripped.Sequence", # PG IDs
                                    expression_col = "Precursor.Quantity",  # PG quantity. Normalized one is also ok。
                                    shrink_ident_cols  = c("Protein.Group","Genes"),  # Gene names
                                    extend_ident_cols =  c("Global.Q.Value"),  # optional, some identification info.
                                    remove_sample_prefix = FALSE,  # remove prefix in sample_col
                                    remove_sample_suffix = FALSE   # remove suffix in sample_col
                                    )

## ----diann_wide2--------------------------------------------------------------
DT::datatable(head(Diann_pep_wided), options = list(scrollX = T,pageLength = 6))

## ----diann_pe-----------------------------------------------------------------
ecols = 2:7

pe <- make_pe_parse(Diann_pep_wided, columns = ecols, 
                    mode = "delim", 
                    remove_suffix = T ## remove suffix
                    )


## ----spe_load2----------------------------------------------------------------
Spe_repfile <- "./A_spikeIn/Spectronaut_output/DIA_MSStats_Report.xls.gz"
Spe_rep <- read.csv(Spe_repfile,sep = "\t") #
dim(Spe_rep) # 318132 rows, each row is a precursor.

## Store both Precursor-, Peptides- and PG-level quality
DT::datatable(head(Spe_rep) %>% mutate_if(is.character, utf8::utf8_encode), options = list(scrollX = T,pageLength = 6))

## ----spe_resha2---------------------------------------------------------------
Spe_rep_wided2 = reshape_long2wide(Spe_rep,
                                   sample_col = "R.FileName",      # the column labeling sample names
                                   feature_col = "PEP.StrippedSequence", # peptide sequence
                                   expression_col = "PEP.Quantity",      # peptide quantity.
                                   shrink_ident_cols  = c("PG.Genes","PG.ProteinAccessions"),  # Gene names and protien IDs
                                   remove_sample_prefix = TRUE,    # remove prefix in sample_col
                                   remove_sample_suffix = FALSE  
                                   )

## ----Spe_pe-------------------------------------------------------------------
ecols = 2:7

pe <- make_pe_parse(Spe_rep_wided2, columns = ecols, 
                    mode = "delim", 
                    remove_suffix = T ## remove suffix
                    )
colData(pe)

## ----sessioninfo, echo=FALSE--------------------------------------------------
sessionInfo()

