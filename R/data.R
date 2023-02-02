#' Silicosis_pg - proteinGroups data of a silicosis mouse model
#'
#' The Silicosis_pg is proteinGroups result of a proteomics study on silica-exposed mouse model
#' reported by Wang, Zhang, Liu, et al.
#'
#' @format A data.frame with 10718 rows and 55 variables:
#' \describe{
#'   \item{Protein.IDs}{Uniprot IDs}
#'   \item{Majority.protein.IDs}{Uniprot IDs of major protein(s) in the protein group}
#'   \item{Peptide.counts..all.}{Total peptides counts of each protein}
#'   \item{Peptide.counts..razor.unique.}{Razor or unique peptides counts of each protein}
#'   \item{Peptide.counts..unique.}{Unique peptides counts of each protein}
#'   \item{Protein.names}{Full protein names}
#'   \item{Gene.names}{Gene name}
#'   \item{Fasta.headers}{Header as present in the Uniprot fasta file}
#'   \item{Peptides}{Number of peptides identified for this protein group}
#'   \item{Razor...unique.peptides}{Number of peptides used for the quantification of this protein group}
#'   \item{Unique.peptides}{Number of peptides identified which are unique for this protein group}
#'   \item{iBAQ columns}{intensity-based absolute-protein-quantification}
#'   \item{LFQ.intensity columns}{LFQ normalized mass spectrometry intensity, A.U.}
#'   \item{Only.identified.by.site}{The protein is only identified by a modification site if marked ('+')}
#'   \item{Reverse}{The protein is identified in the decoy database if marked ('+')}
#'   \item{Potential.contaminant}{The protein is a known contaminant if marked ('+')}
#' }
#' @source Wang M, Zhang Z, Liu J, et al.
#' Gefitinib and fostamatinib target EGFR and SYK to attenuate silicosis:
#' a multi-omics study with drug exploration. Signal Transduct Target Ther. 2022;7(1):157.
#' \href{https://www.nature.com/articles/s41392-022-00959-3}{doi: 10.1038/s41392-022-00959-3}.
#'
"Silicosis_pg"


#' Silicosis_peptide - peptide quantity data of a silicosis mouse model
#'
#' The Silicosis_peptide is peptide quantity of a proteomics study on silica-exposed mouse model
#' reported by Wang, Zhang, Liu, et al.
#'
#' @format A data.frame with 143528 rows and 29 variables:
#' \describe{
#'   \item{Sequence}{Peptide Sequence}
#'   \item{Missed.cleavages}{Missed cleavages number in peptide}
#'   \item{Proteins}{The uniprot ID of proteins that peptide affiliated to}
#'   \item{Leading.razor.protein}{The leading razor protein uniprot ID that peptide affiliated to}
#'   \item{Gene.names}{Gene name of the proteins that peptide affiliated to}
#'   \item{PEP}{Posterior error probability of the identification}
#'   \item{Score}{Andromeda score in Maxquant}
#'   \item{Intensity columns}{Raw intensity of peptides}
#'   \item{Reverse}{The peptide is identified in the decoy database if marked ('+')}
#'   \item{Potential.contaminant}{The peptide from a known contaminant protein if marked ('+')}
#' }
#' @source Wang M, Zhang Z, Liu J, et al.
#' Gefitinib and fostamatinib target EGFR and SYK to attenuate silicosis:
#' a multi-omics study with drug exploration. Signal Transduct Target Ther. 2022;7(1):157.
#' \href{https://www.nature.com/articles/s41392-022-00959-3}{doi: 10.1038/s41392-022-00959-3}.
#'
"Silicosis_peptide"



#' Silicosis_phos - phosphoproteome result of a silicosis mouse model
#'
#' The enriched phosphopeptide quantity data on a silica-exposed mouse model
#' reported by Wang, Zhang, Liu, et al.
#'
#' @format A data.frame with 21249 rows and 40 variables:
#' \describe{
#'   \item{Proteins}{The proteins uniprot IDs that phosphopeptide affiliated to}
#'   \item{Protein}{The first protein uniprot ID}
#'   \item{Gene.names}{Gene names of the proteins that phosphopeptide affiliated to}
#'   \item{Localization.prob}{The localization probability of the phosphorylated site}
#'   \item{PEP}{Posterior error probability of the identification}
#'   \item{Score}{Andromeda score of peptide identification}
#'   \item{Score.for.localization}{Andromeda score of the MS/MS spectrum used for calculating the localization score for this site}
#'   \item{Amino.acid}{The phosphorylated Amino acid}
#'   \item{Intensity columns}{Raw intensity of phosphopeptide}
#'   \item{Reverse}{The peptide is identified in the decoy database if marked ('+')}
#'   \item{Potential.contaminant}{The peptide from a known contaminant protein if marked ('+')}
#'   \item{Positions}{The possible phosphorylation positions in peptides}
#'   \item{Position}{The most possible position with phosphorylation}
#' }
#' @source Wang M, Zhang Z, Liu J, et al.
#' Gefitinib and fostamatinib target EGFR and SYK to attenuate silicosis:
#' a multi-omics study with drug exploration. Signal Transduct Target Ther. 2022;7(1):157.
#' \href{https://www.nature.com/articles/s41392-022-00959-3}{doi: 10.1038/s41392-022-00959-3}.
#'
"Silicosis_phos"


#' Silicosis_ExpDesign - experiment design of silicosis mouse proteome
#'
#' The ExpDesign table of \code{Silicosis_pg} annotating samples information.
#'
#' @format A data.frame with 20 rows and 4 variables:
#' \describe{
#'   \item{label/ID}{The label of experiment groups}
#'   \item{condition }{The experiment condition}
#'   \item{replicate}{Replicate number}
#' }
#'
"Silicosis_ExpDesign"
