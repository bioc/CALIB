\name{getColClasses}
\alias{getColClasses}
\title{Construct colClasses vector for use within read.rg function}

\description{
 Construct a colClasses vector, an argument in read.table used in read.rg
}

\usage{
 getColClasses(cols,...)
}

\arguments{

  \item{cols}{ a character vector of all columns to search against}
  \item{...}{accepts any character lists, character vectors or weight
             functions to match wanted columns against cols}
}

\details{
 This is an internally called function by read.rg to create a colClasses 
 vector used in read.table for fast loading of only required columns in read.rg
}

\value{
 Character. A named vector of classes to assumed for the columns. Possible values are NA
 (when \code{\link{type.convert}} is used), NULL (when the column is skipped). 
}

\references{ \code{getColClasses} in limma package}

\author{ Hui Zhao}

\examples{

allcnames <- c("Block","Column","Row","Name","ID", "F635 Mean","F532 Mean","Flag","Autoflag")
Annotation <- c("Block","Column")
Columns <- list(R="F635 Mean",G="F532 Mean")

getColClasses(allcnames, Annotation, Columns)
}

\keyword{file}
