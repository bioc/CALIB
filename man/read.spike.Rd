\name{read.spike}
\alias{read.spike}

\title{ Read SpikeList from a RGList\_CALIB and Concentration files}

\description{
    Reads a SpikeList from a given RGList\_CALIB object and user specified concentration files
}

\usage{
read.spike(RG, file = NULL, path = NULL, ext = NULL, sep = "\t", conccol,
           regexpcol, different = FALSE, ...)
}

\arguments{
  \item{RG}{ a \code{\link[CALIB:RGList_CALIB-class]{RGList_CALIB}} object.}

  \item{file}{ character string giving the name of the concentration file.}

  \item{path}{ character string giving the directory containing the file.
        Can be omitted if the file is in the current working directory.}

  \item{ext}{ character string giving optional extension to be added to
        each file name.}

  \item{sep}{ field separator character.}

  \item{conccol}{ list with fields \code{RConc}, \code{GConc} giving the column names
       to be used for red and green concentration in the concentration file.}

  \item{regexpcol}{ character vector giving regular expressions in the 
       concentration file.}

  \item{different}{ a logical value. \code{TRUE} means that different arrays use 
        different spikes. So every array should have one concentration file. 
        \code{FALSE} means that different arrays use the same spike and only one 
        concentration file is needed. Default is \code{FALSE}.}

  \item{\dots}{ any other arguments are passed to \code{\link[utils]{read.table}}. }
}

\details{
  This is the function to generate SpikeList in the CALIB package. SpikeList 
  contains all the information for the spikes. This function exacts foreground 
  and background intensities and spot area from RGList\_CALIB, which is generated from 
  funtion \code{\link{read.rg}}. Also this funtion reads concentrations 
  for each spike from a user-specified concentration file (or more than one 
  concentration files if different arrays use different spikes).

  For the concentration file, it should contain the following columns:
  regular expression, red and green concentrations and spike type.
  Spike type should be in the set of \code{Calibration},\cr \code{Ratio}
  and \code{Negative}.

  See the CALIB User's Guide for the example of this function.
}

\value{
  An \code{\link[CALIB:SpikeList-class]{SpikeList}} object containing the components
 
  \item{R}{ matrix containing the red channel foreground intensities for
        each spot for each array.}

  \item{G}{ matrix containing the green channel foreground intensities
        for each spot for each array.}

  \item{Rb}{ matrix containing the red channel background intensities for
        each spike for each array.}

  \item{Gb}{ matrix containing the green channel background intensities
        for each spike for each array.}

  \item{RArea}{ matrix containing the red spot area for each spike for each array.}

  \item{GArea}{ matrix containing the green spot area for each spike for each array.}

  \item{RConc}{ matrix containing the red concentration for each spike for each array.}

  \item{GConc}{ matrix containing the green concentration for each spike for each array.}

  \item{genes}{ data frame containing annotation information about the probes, 
        for example, spike types and IDs and spatial positions on the array}
}

\references{ ~put references to the literature/web site here ~ }

\author{ Hui Zhao}

\keyword{ file }


