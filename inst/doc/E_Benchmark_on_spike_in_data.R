## ----workplace,eval=FALSE-----------------------------------------------------
#  # Download OmicsExample and uncompress firstly
#  # Change the path to Omics Example
#  example_path = "the/path/to/OmicsExample"
#  knitr::opts_knit$set(
#     root.dir = example_path
#   )
#  knitr::opts_chunk$set(
#     collapse = TRUE,
#     comment = "#>",
#     message = FALSE,
#     warning = FALSE
#  )
#  
#  library(DEP2)
#  library(magrittr)
#  library(dplyr)
#  library(ggVennDiagram)
#  library(patchwork)
#  library(ggplot2)

## ----workplace0,include=FALSE,eval=TRUE---------------------------------------
# Download Example data and set workplace
example_path = "../example_data/OmicsExample/"
knitr::opts_knit$set(
   root.dir = example_path
 )
knitr::opts_chunk$set(
   collapse = TRUE,
   comment = "#>",
   message = FALSE,
   warning = FALSE
)

# load required package
library(DEP2)
library(magrittr)
library(dplyr)
library(ggVennDiagram)
library(patchwork)
library(ggplot2)

## ----pg_dep2,results='hide',message=FALSE,warning=FALSE-----------------------
mq_pgfile <- "./A_spikeIn/MaxQuant_output/proteinGroups.txt.gz"
pg <- read.csv(mq_pgfile,sep = "\t")

## Generate unique names and ids
unique_pg <- make_unique(pg, 
                         name = "Gene.names", #gene 
                         ids = "Protein.IDs"   #protein
                         )

## Convert se with expdesign table
ecols <- grep("LFQ.intensity.", colnames(unique_pg))
se <- make_se_parse(unique_pg, columns = ecols, mode = "delim")

## Filter & normalize & Impute
filt <- filter_se(se,thr = 0, fraction = 0.5, 
                  filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" & Peptides > 1)
norm <- normalize_vsn(filt)
set.seed(10000)
imp <- DEP2::impute(norm,fun = "RF")

## Differential test
diff <- test_diff(imp,type = "control", control = "A1", fdr.type = "BH")
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

## ----pep_dep,results='hide',message=FALSE,warning=FALSE-----------------------
mq_pepfile <- "./A_spikeIn/MaxQuant_output/peptides.txt.gz"
mq_pep <- read.csv(mq_pepfile,sep = "\t")

## Construct SummarizedExperiement
ecols <- grep("Intensity.", colnames(mq_pep), value = TRUE)
pe = make_pe_parse(mq_pep, columns = ecols,   # columns is the 'Intensity' colunmns
                   mode = "delim", sep = "_", remove_prefix = TRUE)

##  Filter & Impute & Normalize
pe = filter_pe(pe, thr = 0, fraction = 0.6, 
               filter_formula = ~ Reverse != '+' & Potential.contaminant !="+" )
pe <- normalize_pe(pe,method = "vsn", i = "peptideRaw")
set.seed(10000)
pe <- impute_pe(pe, fun = "bpca", name = "peptideImp", i = "peptideNorm")

## Aggregation
pe <- aggregate_pe(pe, aggrefun = "RobustSummary", fcol = "Proteins", reserve = "Gene.names", i = "peptideImp")
se_pep <- pe2se(pe, names = "Gene.names", ids = "smallestProteingroups")
filt_pep <- filter_se(se_pep,filter_formula = ~.n > 1) 
norm_pep <- normalize_vsn(filt_pep)


## Differential test
diff_pep <- test_diff(norm_pep,type = "control", control = "A1", fdr.type = "BH")
dep_pep <- add_rejections(diff_pep, alpha = 0.01)

## ----pg_perseus,results='hide',message=FALSE,warning=FALSE--------------------
## The differetial test result of Perseus
Perseus_resfile <- "./A_spikeIn/MaxQuant_output/Perseus_test_res.txt.gz"
con <- file(Perseus_resfile, "r")
Perseus_resfile = readLines(con)
close(con)

## Remove the lines recording variable types.
Perseus_resfile = Perseus_resfile[-(2:3)]

## Read data
Perseus_res <- read.csv(text = Perseus_resfile,sep = "\t")

## Clean data
Perseus_res$p.value = 10^-(Perseus_res$X.Log.Student.s.T.test.p.value.B1_A1)
Perseus_res$q.value = Perseus_res$Student.s.T.test.q.value.B1_A1
Perseus_res$L2FC = Perseus_res$Student.s.T.test.Difference.B1_A1
Perseus_res$t.statistic = Perseus_res$Student.s.T.test.Test.statistic.B1_A1
Perseus_res <- make_unique(Perseus_res, "Gene.names", "Protein.IDs", delim = ";")
Perseus_res <- dplyr::filter(Perseus_res,Peptides > 1)

## Convert a SE object for subsequent comparison
ecols <- 1:6
Perseus_res2 = Perseus_res
Perseus_res2[ecols] <- 2^Perseus_res2[ecols]
dep_Perseus = make_se_parse(Perseus_res2,columns = 1:6,mode = "delim")

## ----dia_dep,results='hide',message=FALSE,warning=FALSE-----------------------
Spe_repfile <- "./A_spikeIn/Spectronaut_output/DIA4_MSStates_Report_new_noimpu.xls"
Spe_rep <- read.csv(Spe_repfile,sep = "\t") 

## Reshape long table 
Spe_rep_wided <- reshape_long2wide(Spe_rep,
                                  sample_col = "R.FileName",      # the column labeling sample names
                                  feature_col = "PG.ProteinGroups", # PG IDs
                                  expression_col = "PG.Quantity",   # PG quantity.
                                  shrink_ident_cols  = c("PG.Genes","PG.Qvalue","PG.RunEvidenceCount"),  # Gene names
                                  remove_sample_prefix = TRUE,    # remove prefix in sample_col
                                  remove_sample_suffix = FALSE  
                                  )

## Unique names
unique_spe_pg <- make_unique(Spe_rep_wided, 
                             name = "PG.Genes",  #gene 
                             ids = "PG.ProteinGroups"   #protein
                             )

## Convert SE with expdesign table
ecols <- 2:7 # Expression col in the Spectronaut output
se_dia <- make_se_parse(unique_spe_pg, 
                    columns = ecols, 
                    mode = "delim",
                    remove_prefix = TRUE,   ## remove the file prefix 
                    remove_suffix = FALSE   ## no suffix in this case
                    )

## Filter. Keep the protein which have at least two Evidence in one sample
Keep_proteins <- Spe_rep_wided$PG.RunEvidenceCount %>% 
  sapply(function(x){
    max(strsplit(x,split = ";")[[1]]) >1
  }) 
filt_dia <- se_dia[Keep_proteins,]

## Normalize & impute
norm_dia <- normalize_vsn(filt_dia)
imp_dia <- DEP2::impute(filt_dia, "bpca")

## Differential test
diff_dia <- test_diff(imp_dia, type = "control", control = "A5", fdr.type = "Storey's qvalue")
dep_dia <- add_rejections(diff_dia, alpha = 0.05, lfc = 1)

## -----------------------------------------------------------------------------
(plot_coverage(se) + ggplot2::labs(title = "DDA")) + (plot_coverage(se_dia)  + ggplot2::labs(title = "DIA")) +
  (plot_coverage(filt) + ggplot2::labs(title = "DDA-filtered")) + (plot_coverage(filt_dia)  + ggplot2::labs(title = "DIA-filtered")) 

## ----message=FALSE,warning=FALSE----------------------------------------------
features_pg = rowData(dep)$ID
features_pep = rowData(dep_pep)$ID
features_perseus = Perseus_res$ID
features_dia = rowData(dep_dia)$ID

ggVennDiagram(list("DDA-Protein-DEP2" = features_pg,
                   "DDA-Protein-Perseus" = features_perseus,
                   "DDA-Peptide-DEP2" = features_pep,
                   "DIA-Protein-DEP2" = features_dia
                   ),
              set_size = 2.5,
              label_size = 3)

## ----message=FALSE,warning=FALSE----------------------------------------------
protein_info = read.csv("./A_spikeIn/protein_info.csv")
protein_info2 = protein_info %>% dplyr::select(ID,spe,ox)

PG_df <- get_df_wide(dep) %>% as.data.frame %>%  
  mutate( meanQuant = rowMeans(dplyr::select(., matches("^(A1|B1)_(1|2|3)$"))) ) %>% 
  dplyr::select(name, ID, L2FC = B1_vs_A1_diff, pval = B1_vs_A1_p.val, qval = B1_vs_A1_p.adj, t.statistic = B1_vs_A1_t.stastic, meanQuant) %>% 
  left_join(protein_info2) %>% 
  cbind(assay(dep))

PEP_df <- get_df_wide(dep_pep) %>% as.data.frame %>%  
  mutate( meanQuant = rowMeans(dplyr::select(., matches("^(A1|B1)_(1|2|3)$"))) ) %>% 
  dplyr::select(name, ID, 
                L2FC = B1_vs_A1_diff, pval = B1_vs_A1_p.val, qval = B1_vs_A1_p.adj, t.statistic = B1_vs_A1_t.stastic, meanQuant) %>% 
  left_join(protein_info2) %>% 
  cbind(assay(dep_pep))

perseus_df <- Perseus_res %>%  
  mutate( meanQuant = rowMeans(dplyr::select(., matches("^(LFQ.).*_(1|2|3)$"))) ) %>% 
  dplyr::select(name, ID, L2FC, pval = p.value, qval = q.value, t.statistic, meanQuant) %>%
  left_join(protein_info2) %>% 
  cbind(dplyr::select(Perseus_res,starts_with("LFQ.")) )

DIA_df <- get_df_wide(dep_dia) %>% as.data.frame  %>% 
  mutate( meanQuant = rowMeans(dplyr::select(., matches("^(A5|B5)_(1|2|3)$"))) ) %>% 
  dplyr::select(name, ID, L2FC = B5_vs_A5_diff, pval = B5_vs_A5_p.val, qval = B5_vs_A5_p.adj, t.statistic = B5_vs_A5_t.stastic, meanQuant) %>% 
  left_join(protein_info2) %>% 
  cbind(assay(dep_dia)) 

df_list <- list("DDA-Protein-DEP2" = PG_df, 
                "DDA-Peptide-DEP2" = PEP_df, 
                "DDA-Protein-Perseus" = perseus_df, 
                "DIA-DEP2" = DIA_df)

cols = c("ID", "spe", "L2FC", "meanQuant","pval", "qval")
res_table = df_list %>% 
  lapply(., function(x){
    dplyr::select(x,!!!syms(cols)) %>% distinct(., ID, .keep_all = T)
  }) %>% 
  Map(function(x, method){
    x$method <- method
    x
  },., names(.)) %>% 
  do.call(rbind, .) 

intersect_id = df_list %>% lapply(function(x) x$ID) %>% Reduce(intersect, .)


## Plot L2FC boxplot. Dash line is the actual fold change
p_lfc <- ggplot(res_table, aes( y= L2FC, x = spe, fill = method)) +
  geom_boxplot(
    color = "black", size = 0.5, 
    outlier.color = "black", outlier.alpha = 0.3, outlier.size = 0.5, outlier.shape = "."
    ) +
  ylab("Log2 fold change") +
  geom_segment(y= log2(3),yend= log2(3), x = 0.5, xend= 1.5, 
               color = "darkred", linetype = "dashed", lwd = 0.8, alpha = 0.6) +
  geom_segment(y= -log2(2),yend= -log2(2), x = 2.5, xend= 3.5, 
               color = "darkred", linetype = "dashed", lwd = 0.8, alpha = 0.6) +
  geom_segment(y= 0,yend= 0, x = 1.5,xend= 2.5, 
               color = "darkred", linetype = "dashed", lwd = 0.8, alpha = 0.6) +
  scale_y_continuous(limits = c(-3, 4)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_blank(),
        # axis.ticks.length = unit(0, "pt"),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white")
        )
p_lfc

## ----message=FALSE,warning=FALSE----------------------------------------------
res_table2 = dplyr::filter(res_table, ID %in% intersect_id)
res_table2$L2FC = ifelse(res_table2$spe == "HUMAN", 
                          res_table2$L2FC ,
                          ifelse(res_table2$spe == "YEAST", 
                                 res_table2$L2FC + log2(2), 
                                 res_table2$L2FC - log2(3)))

p_lfcerr <- ggplot(res_table2, aes( y= L2FC, x = method, fill = spe)) +
  geom_boxplot(
    color = "black", size = 0.5, 
    outlier.color = "black", outlier.alpha = 0.3, outlier.size = 0.5, outlier.shape = "."
    ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_fill_manual(values = c(HUMAN = "#41b349", YEAST = "#2775b6", ECOLI = "#a6522c")) + 
  scale_y_continuous(name = "Fold change error(log2)", limits = c(-1.5, 1.5)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_blank(),
        # axis.ticks.length = unit(0, "pt"),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white")
        )
p_lfcerr

## ----message=FALSE,warning=FALSE----------------------------------------------
res_table3 = dplyr::filter(res_table, ID %in% intersect_id)
res_table3$L2FC[res_table3$L2FC < -3.5] = -3.5
res_table3$L2FC[res_table3$L2FC > 3.5] = 3.5

p_lfctrend <- ggplot(res_table3, aes(x = meanQuant, y = L2FC, color = spe)) + 
  geom_point(shape = 20, alpha = 0.2, size = 0.8) +
  geom_line(stat="smooth",method = "loess", formula = y ~x, span=0.3, size = 0.7, alpha = 0.5) +
  scale_color_manual(values = c(HUMAN = "#41b349", YEAST = "#2775b6", ECOLI = "#a6522c")) +
  scale_y_continuous(limits = c(-3.5,3.5), name = "Log2 fold change") +
  xlab("Mean abundance") +
  facet_wrap(~method,scales = 'free') +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", colour = NA))
p_lfctrend

## ----message=FALSE,warning=FALSE----------------------------------------------
dep_list = list(DDA_Protein_DEP2 = dep,
                DDA_Peptide_DEP2 = dep_pep,
                DDA_Protein_Perseus = dep_Perseus,
                DIA_DPE2 = dep_dia)

CV_df = dep_list %>% 
  lapply(function(x){
    assa = 2^assay(x) %>% as.data.frame()
    cd = colData(x)
    
    coef_variation <- function(x) {
      coef = sd(x)/mean(x)
    }
    
    CVs = assa %>% tibble::rownames_to_column() %>% 
      cbind(., id = rowData(x)$ID) %>%
      tidyr::gather("label", "Intensity", -c(rowname,id)) %>%
      dplyr::left_join(., data.frame(cd), by = "label") %>%
      dplyr::group_by(id, condition) %>% dplyr::summarise(cvs = coef_variation(Intensity)) %>%
      dplyr::group_by(condition) %>% dplyr::mutate(condition_median = median(cvs)) %>%
      dplyr::ungroup()
    
  }) %>%
  Map(function(x, method){
    x$method <- method
    x
  },., names(.))  %>%
  do.call(rbind,.)

CV_df2 = dplyr::filter(CV_df, id %in% intersect_id) %>% 
  arrange(id) %>%
  mutate(condition = gsub("[0-9]","",condition))


P_A <- wilcox.test(
  dplyr::filter(CV_df2, method == "DDA_Protein_DEP2", condition == "A") %>% .$cvs, 
  dplyr::filter(CV_df2, method == "DDA_Protein_Perseus", condition == "A") %>% .$cvs,
  paired = T)

P_B <- wilcox.test(
  dplyr::filter(CV_df2, method == "DDA_Protein_DEP2", condition == "B") %>% .$cvs, 
  dplyr::filter(CV_df2, method == "DDA_Protein_Perseus", condition == "B") %>% .$cvs,
  paired = T)

library(ggsignif)
p_cvs <- ggplot(CV_df2, aes(x = condition, y = cvs, fill = method)) +
    geom_boxplot(color = "black", 
                 outlier.size = 0.2,outlier.alpha = 0.2, position = "dodge2") +
    geom_signif(annotation = ifelse(P_A$p.value < 0.001, "p<0.001", P_A$p.value),
                y_position = 1.6, xmin = 0.91, xmax = 1.09, 
                tip_length = c(0.3, 0.03)
    ) +
    geom_signif(annotation = ifelse(P_B$p.value < 0.001, "p<0.001", P_B$p.value),
                y_position = 1.65, xmin = 1.91, xmax = 2.09, 
                tip_length = c(0.45, 0.02)
    ) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_blank(),
          # axis.ticks.length = unit(0, "pt"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white")
    )

p_cvs

## ----message=FALSE,warning=FALSE----------------------------------------------
library(pROC)
cols = c("ID", "spe", "L2FC", "pval", "qval", "t.statistic")
df_list %>% 
  lapply(function(x){
    dplyr::select(x,!!!syms(cols)) %>% 
      distinct(., ID, .keep_all = T) %>% 
      dplyr::filter(ID %in% intersect_id)
  }) %>% lapply(function(x){
    pROC::roc(ifelse(x$spe == "HUMAN", F, T), abs(x$qval))
  }) -> roc_res

roc_res %>% lapply(function(x){
  df <- data.frame(specificities = x$specificities, sensitivities = x$sensitivities, qval = x$thresholds)
}) %>%
  Map(function(x, method){
    x$method <- method
    x
  },., names(.)) %>% 
  do.call(rbind,.) -> roc_df

roc_res %>% sapply(function(x){
  as.numeric(x$auc)
}) -> aucs 
aucs <- paste(names(aucs),round(aucs,3),sep = "  ", collapse = ", ")

ggplot(roc_df) +
  geom_line( aes(x=1 - specificities, y=sensitivities, color=method)) + 
  coord_equal() +
  xlab("") + ylab("") +
  scale_x_continuous(limits = c(0, 0.2), n.breaks = 3, labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(limits = c(0.8, 1), n.breaks = 3, labels = function(x) sprintf("%.2f", x)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") -> p_zoom

ggplot(roc_df) +
  geom_line( aes(x=1 - specificities, y=sensitivities, color=method)) + 
  xlim(c(0, 1)) + ylim(c(0, 1)) +
  coord_equal() +
  xlab("1 - Specificity") + ylab("Sensitivity") +
  ggtitle("ROC Curves for Colon Cancer Markers") +
  annotate("rect", xmin = 0, xmax = 0.2, ymin = 0.8, ymax = 1, fill = "transparent", color = "grey",alpha = 0.2) +
  annotation_custom(ggplotGrob(p_zoom), xmin = 0.25, xmax = 0.9, ymin = -0.2, ymax = 0.9) +
  theme_bw() +
  theme(panel.grid = element_blank())

## ----message=FALSE,warning=FALSE----------------------------------------------
ggplot(roc_df) +
  geom_point( aes(x=qval, y=1-specificities, color=method), size = 0.5) +
  xlim(c(0, 0.05)) + ylim(c(0, 0.05)) +
  geom_abline(intercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank())

## ----message=FALSE,warning=FALSE----------------------------------------------
res_table$trend = apply(res_table, 1, 
      function(x){
        if( as.numeric(x["qval"]) > 0.05)
          return("Not sig")
        if( as.numeric(x["qval"]) <= 0.001 & as.numeric(x["L2FC"]) > 1.2)
          return("Up")
        if( as.numeric(x["qval"]) <= 0.05 & as.numeric(x["L2FC"]) < -0.6)
          return("Down")
        
        return("Not sig")
      },simplify = F) %>% unlist

res_table$trend = factor(res_table$trend,c("Not sig","Up","Down"))
res_table$spe = factor(res_table$spe,c("ECOLI", "YEAST", "HUMAN"))

ggplot(res_table) + geom_bar(aes(x = trend, fill = spe,y=1),stat="identity") +
  scale_fill_manual(values = c(HUMAN = "#41b349", YEAST = "#2775b6", ECOLI = "#a6522c")) +
  facet_grid(~method,scales = 'free') + 
  ylab("Proteingroups") + 
  xlab("Significant") + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = NA)
  )

## ----message=FALSE,warning=FALSE----------------------------------------------
ggplot(res_table) + geom_bar(aes(x = trend, fill = spe,y=1), stat="identity", position = "fill" ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(HUMAN = "#41b349", YEAST = "#2775b6", ECOLI = "#a6522c")) +
  facet_grid(~method,scales = 'free') + 
  ylab("Proteingroups") + 
  xlab("Significant") + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_blank(),
        # axis.ticks.length = unit(0, "pt"),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = NA)
  )

## ----include=FALSE,eval=FALSE-------------------------------------------------
#  library(Biostrings)
#  Biostrings::readAAStringSet("./A_spikeIn/YeastE12_S.cerevisiae_homo_20211201_uniprot-reviewed_yes+taxonomy_9605.fasta") -> fa
#  names = names(fa)
#  protein_info = names(fa) %>% lapply(function(x){
#    ID = stringr::str_extract(x,"(?<=\\|)[^ ]+(?=\\|)")
#    name <- stringr::str_extract(x, "(?<=GN=)[^ ]+")
#    spe <- stringr::str_extract(x, "(?<=_)[^ ]+")
#    ox <- stringr::str_extract(x, "(?<=OX=)[^ ]+")
#  
#    data.frame(ID, name, spe, ox, header = x)
#  })
#  
#  protein_info <- do.call(rbind, protein_info)
#  
#  write.csv(protein_info,"./A_spikeIn/protein_info.csv",row.names = F)

## ----sessioninfo, echo=FALSE--------------------------------------------------
sessionInfo()

