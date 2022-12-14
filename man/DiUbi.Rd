% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{DiUbi}
\alias{DiUbi}
\title{DiUbi - Ubiquitin interactors for different diubiquitin-linkages (UbIA-MS dataset)}
\format{
A data.frame with 4071 observations and 102 variables:
\describe{
  \item{Protein.IDs}{Uniprot IDs}
  \item{Majority.protein.IDs}{Uniprot IDs of major protein(s) in the protein group}
  \item{Protein.names}{Full protein names}
  \item{Gene.names}{Gene name}
  \item{Fasta.headers}{Header as present in the Uniprot fasta file}
  \item{Peptides}{Number of peptides identified for this protein group}
  \item{Razor...unique.peptides}{Number of peptides used for the quantification of this protein group}
  \item{Unique.peptides}{Number of peptides identified which are unique for this protein group}
  \item{Intensity columns (30)}{Raw mass spectrometry intensity, A.U.}
  \item{iBAQ columns (30)}{iBAQ normalized mass spectrometry intensity, A.U.}
  \item{LFQ.intensity columns (30)}{LFQ normalized mass spectrometry intensity, A.U.}
  \item{Only.identified.by.site}{The protein is only identified by a modification site if marked ('+')}
  \item{Reverse}{The protein is identified in the decoy database if marked ('+')}
  \item{Potential.contaminant}{The protein is a known contaminant if marked ('+')}
  \item{id}{The protein group ID}
}
}
\source{
Zhang, Smits, van Tilburg, et al (2017).
An interaction landscape of ubiquitin signaling.
Molecular Cell 65(5): 941-955.
\href{http://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30004-7}{doi: 10.1016/j.molcel.2017.01.004}.
}
\usage{
DiUbi
}
\value{
A data.frame.
}
\description{
The DiUbi dataset contains label free quantification (LFQ)
and intensity-based absolute quantification (iBAQ) data for
ubiquitin interactors of different diubiquitin-linkages,
generated by Zhang et al 2017.
The dataset contains the proteingroups output file from
\href{http://www.maxquant.org}{MaxQuant}.
}
\keyword{datasets}
