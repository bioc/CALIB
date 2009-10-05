\name{read.rg}
\alias{read.rg}
\title{Read RGList\_CALIB from Image Analysis Output Files}
\description{
 Reads an RGList\_CALIB from a series of microarray image analysis output files
}

\usage{
read.rg(files = NULL, source = "generic", path = NULL, ext = NULL, 
        names = NULL, columns = NULL, other.columns = NULL,
        annotation = NULL,wt.fun = NULL, verbose = TRUE, sep = "\t",
        quote = NULL, DEBUG = FALSE, ...)
}

\arguments{
  \item{files}{ character vector giving the names of the files containing
         image analysis output or, for Imagene data, a character
         matrix of names of files. If omitted, then all files with
         extension \code{ext} in the specified directory will be read in
         alphabetical order.}

  \item{source}{ character string specifying the image analysis program which
         produced the output files.  Choices are \code{"generic"},
         \code{"agilent"}, \code{"arrayvision"}, \code{"bluefuse"}, \code{"genepix"},
         \code{"genepix.custom"}, \code{"genepix.median"},\cr \code{"imagene"},
         \code{"quantarray"}, \code{"scanarrayexpress"}, \code{"smd.old"},
         \code{"smd"},\cr \code{"spot"} or \code{"spot.close.open"}.}

  \item{path}{ character string giving the directory containing the files.
         The default is the current working directory.}

  \item{ext}{ character string giving optional extension to be added to
         each file name}

  \item{names}{ character vector of names to be associated with each array as
         column name. Defaults to \code{removeExt(files)}.}

  \item{columns}{ list with fields \code{R}, \code{G}, \code{Rb}, \code{Gb},
        \code{RArea} and \code{GArea} giving the column names to be used for
        red foreground, green foreground, red background, green background,
        red area and green area respectively. Or, in the case of Imagene data,
        a list with fields \code{f} and \code{b}. This argument is optional
        if \code{source} is specified, otherwise it is required.}

  \item{other.columns}{ character vector of names of other columns to be read
         containing spot-specific information}

  \item{annotation}{ character vector of names of columns containing annotation
         information about the probes}

  \item{wt.fun}{ function to calculate spot quality weights}

  \item{verbose}{ logical, \code{TRUE} to report each time a file is read }

  \item{sep}{ the field separator character }

  \item{quote}{ character string of characters to be treated as quote marks }

  \item{DEBUG}{ a logical value, if \code{TRUE}, a series of echo statements 
         will be printed for each file. Details on the file, skip, 
         and selected columns in a colClasses format for read.table 
         will be displayed.}

  \item{\dots}{ any other arguments are passed to
    \code{\link[utils]{read.table}}.} 

}

\details{
  This is the main data input function for CALIB package. It has the similar
  usage as the read.maimages function in limma package. The output of the 
  function is an \code{RGList_CALIB} object. However,there
  are two more fields - \code{$RArea} and \code{$GArea} than \code{RGList}
  object in limma package. These two fields contain spot area of each color.
  More details see \code{read.maimages} in limma package.
}

\value{
  An \code{\link[CALIB:RGList_CALIB-class]{RGList_CALIB}} object
  containing the components 
 
  \item{R}{ matrix containing the red channel foreground intensities for
        each spot for each array.}

  \item{G}{ matrix containing the green channel foreground intensities
        for each spot for each array.}

  \item{Rb}{ matrix containing the red channel background intensities for
        each spot for each array.}

  \item{Gb}{ matrix containing the green channel background intensities
        for each spot for each array.}

  \item{RArea}{ matrix containing the red spot area for each spot for each array.}

  \item{GArea}{ matrix containing the green spot area for each spot for each array.}

  \item{weights}{ spot quality weights, if \code{wt.fun} is given}

  \item{other}{ list containing matrices corresponding to \code{other.columns} if given}

  \item{genes}{ data frame containing annotation information about the
        probes, for example gene names and IDs and spatial positions
        on the array, currently set only if \code{source} is \code{"agilent"},
        \code{"genepix"} or \code{source="imagene"} or if the \code{annotation}
        argument is set}

  \item{targets}{ data frame with column \code{FileName} giving the names of the
        files read}

  \item{source}{ character string giving the image analysis program name}

  \item{printer}{ list of class \code{PrintLayout}, currently set only if
        \code{source="imagene"}}
}

\references{ \code{\link[limma]{read.maimages}} in limma package}

\author{ Hui Zhao}

\seealso{ 'read.rg' is based on \code{\link[utils]{read.table}} in the base package}

\examples{

 #  Read all .gpr files from current working directory.

 # files <- dir(pattern="*\\.gpr$")
 # RG <- read.rg(files,"genepix")
}

\keyword{ file }% at least one, from doc/KEYWORDS





