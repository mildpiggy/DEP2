# data(peptide3)

Silicosis_peptide <- read.csv("example/silicosis_peptides.txt",sep = "\t")
head(Silicosis_peptide)
colnames(Silicosis_peptide)
retain_cols <- c(which(colnames(Silicosis_peptide) %in% c("Sequence", "Length", "Missed.cleavages", "Mass", "Proteins", "Leading.razor.protein",
                                                          "Gene.names", "Protein.names", "Charges", "PEP", "Score")),
                 99:119,120:122)
Silicosis_peptide = Silicosis_peptide[,retain_cols]
usethis::use_data(Silicosis_peptide, overwrite = T, compress = "xz")


Silicosis_pg <- read.csv("example/silicosis_pGs.txt",sep = "\t")
colnames(Silicosis_pg)
Silicosis_pg <- Silicosis_pg[,c(1:12,156:195,217:219)]
usethis::use_data(Silicosis_pg, overwrite = T,compress = "xz")

Silicosis_phos <- read.csv("example/sillicosis_Phospho (STY)Sites.txt",sep = "\t")
head(Silicosis_phos)
colnames(Silicosis_phos)
retain_cols <- c(which(colnames(Silicosis_phos) %in% c("Proteins", "Positions.within.proteins", "Leading.proteins", "Protein", "Protein.names", "Gene.names",
                                                          "Gene.names", "Fasta.headers", "Localization.prob", "PEP", "Score", "Score.for.localization", "Amino.acid",
                                                          "Charge", "Intensity")),
                 130:149,
                 which(colnames(Silicosis_phos) %in% c("Reverse","Potential.contaminant", "id", "Protein.group.IDs", "Positions", "Position")))
Silicosis_phos <- Silicosis_phos[,retain_cols]
usethis::use_data(Silicosis_phos, overwrite = T, compress = "xz")

# to reduce store space of .rds
tools::resaveRdaFiles("data/Silicosis_phos.rda")
tools::checkRdaFiles("data/")


# colnames(Silicosis_pg)
# unique_pg <- make_unique(Silicosis_pg, names = "Gene.names", ids = "Protein.IDs")
# ecols <- grep("LFQ.intensity", colnames(unique))
# se <- make_se_parse(unique, columns = ecols, mode = "delim", sep = "_")
# plot_frequency(se)
# DEP2::plot_coverage(se)
# filtered <- filter_se(se, thr = 0, filter_formula = ~ Reverse != "+" & Potential.contaminant!="+")
# plot_frequency(filtered)
# DEP2::plot_coverage(filtered)






