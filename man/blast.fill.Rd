% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blast.fill}
\alias{blast.fill}
\title{Local BLAST results with added content.}
\format{
A data frame with 72 rows and 19 variables:
\describe{
  \item{probe}{probe sequence}
  \item{probe.length}{probe sequence's length}
  \item{Qid}{query identification number}
  \item{Qstart}{query start position}
  \item{Qend}{query end position}
  \item{Rgi}{subject GenInfo Identifier number}
  \item{Racc}{subject NCBI accession number}
  \item{Rtitle}{subject title}
  \item{Rtaxid}{subject taxon identificator}
  \item{Rstart}{subject start position}
  \item{Rend}{subject end position}
  \item{alig.length}{length of alignment}
  \item{mismatch}{amount of mismatches}
  \item{gaps}{amount of gaps}
  \item{ident.number}{amount of identical positions}
  \item{score}{alignment score}
  \item{bitscore}{alignment bitscore}
  \item{Evalue}{alignment e-value}
  \item{Qcover}{query coverage, \%}
  }
}
\usage{
blast.fill
}
\description{
Result of BLAST of 5 probes against local database of target nucleotide sequences
of Chlamydia pneumoniae. Local BLAST was performed with blast_local () function.
Subjects' Genbank Identifiers are added with fill_blast_result () function.
}
\keyword{datasets}
