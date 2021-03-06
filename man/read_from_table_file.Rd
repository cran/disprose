% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common.R
\name{read_from_table_file}
\alias{read_from_table_file}
\title{Read table file}
\usage{
read_from_table_file(
  file,
  choose.columns = FALSE,
  column.names,
  select = FALSE,
  select.column.name,
  select.val,
  unique = FALSE,
  sep = ";",
  header = TRUE
)
}
\arguments{
\item{file}{character; file name and path}

\item{choose.columns}{logical; return chosen columns only}

\item{column.names}{character; vector of name of columns that are chosen to be returned}

\item{select}{logical; return only rows that contain selected values in one column}

\item{select.column.name}{character; name of column that contains selected values}

\item{select.val}{vector of values that define rows that should be returned}

\item{unique}{logical; delete duplicated rows}

\item{sep}{character; the field separator character}

\item{header}{logical; files contain the names of the variables as its first line}
}
\value{
Data frame with file content, optionally trimmed.
}
\description{
Read table file and selects the required rows and columns.
}
\details{
This function reads table files and returns data frame with selected
rows (only rows with specified values) and columns.
Also duplicated rows may be deleted.

\code{column.names} and  \code{select.column.name} must be exact column names as in data frame.
}
\examples{
mydata <- data.frame (N = 1:10, letter = c(rep ("A", 5), rep ("B", 4), "C"),
                   num = c(1, rep(1:4, 2), 5))
t1<-tempfile()
write.table (mydata, t1, sep = ";")
read_from_table_file (file = t1)
read_from_table_file (file = t1, select = TRUE, select.column.name = "letter",
                     select.val = c("A", "C"))
read_from_table_file (file = t1, select = TRUE, select.column.name = "letter",
                     select.val = c("A", "C"), unique=TRUE, choose.columns = TRUE,
                     column.names = c("letter", "num"))
read_from_table_file (file = t1, select = TRUE, select.column.name = "letter",
                     select.val = c("A", "C"), unique = TRUE, choose.columns = TRUE,
                     column.names = c("N", "num"))
read_from_table_file (file = t1, select = TRUE, select.column.name = "letter",
                     select.val = c("A", "C"), unique = TRUE, choose.columns = TRUE,
                     column.names = c("letter", "N"))
file.remove (t1)

}
\author{
Elena N. Filatova
}
