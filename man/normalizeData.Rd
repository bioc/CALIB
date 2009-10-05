\name{normalizeData}
\alias{normalizeData}
\title{ Normalization: estimation of absolute expression levels}

\description{
   estimates absolute expression levels for each combination of a gene 
   and a tested biological condition.
}

\usage{
normalizeData(RG, parameter, array = array, condition = condition, dye = dye,
              cloneid = cloneid, idcol = idcol)
}

\arguments{
  \item{RG}{ an \code{\link[CALIB:RGList_CALIB-class]{RGList_CALIB}} object}
 
  \item{parameter}{ a
    \code{\link[CALIB:ParameterList-class]{ParameterList}} object} 

  \item{array}{ integer vector specifying the index of the arrays. 
        Has length equal to two times of the number of arrays. }

  \item{condition}{ integer vector specifying the index of the conditions. 
        Has length equal to two times of the number of arrays. }

  \item{dye}{ integer vector specifying the index of the dyes. Has length 
        equal to two times of the number of arrays. }

  \item{cloneid}{ string vector specifying the clone ids of the clones to 
        be normalized. If missing, normalize all the clones.}

  \item{idcol}{ string specifying the column name of clone ids in the 
        genes field of RG. }
 }

\details{
  This function estimates absolute expression levels for each combination of 
  a gene and a tested biological condition from the measured intensity. 
  It accepts measured intensities from \code{RG}.

  The argument \code{parameter} is an object of \code{ParameterList}. The function accepts 
  model parameters from this argument.

  By using this function, for each combination of a gene and a tested biological
  condition, a single absolute expression level fo target is estimated. Therefore,
  specifying the design of experiment is necessary. Namely, the design of array, 
  condition and dye is needed. The three arguments \code{array},\code{condition} 
  and \code{dye} are three numeric vector to indicate the design of array, 
  condition and dye respectively. How to specify these three arguments refer to 
  the example below.

  The function is able to not only estimate all the genes on the slides but also
  estimate any gene on the slides seperately. The argument \code{cloneid} accepts the
  clone ids of which the genes are interested by the user. If \code{cloneid} argument
  is missing, the function will estimate all the genes on the slides.In order to 
  match clone id in the RG, column name which indicates clone ids in RG\$genes
  should be specified by argument \code{idcol}. 
}

\note{
  The main calculation part in this function is done by c++ code.
}

\section{Warning}{
  The function doesn't allow missing clone id. So please check before run the function.
}

\value{
  a numeric matix containing the absolute expression levels. Columns indicate
  different conditions and rows indicate different genes.
}

\references{ Engelen, K., Naudts, B., DeMoor, B., Marchal, K. (2006) 
   A calibration method for estimating absolute expression levels from 
   microarray data. Bioinformatics  22: 1251-1258.}

\author{ Hui Zhao}

\examples{

# load data: RG and parameter
data(RG)
data(parameter)

# define design matrix: two arrays, two condition and color-flip design
array <- c(1,1,2,2)
condition <- c(1,2,2,1)
dye <- c(1,2,1,2)

# specify clone-id column
idcol <- "CLONE_ID"

#data <- normalizeData(RG,parameter,array=array,condition=condition,dye=dye,idcol=idcol)

## only normalize a group of genes
cloneid_interested <- c("250001", "250002", "250003", "250004", "250005")
data <- normalizeData(RG,parameter,array=array,condition=condition,dye=dye,cloneid=cloneid_interested,idcol=idcol)

}

\keyword{ optimize }
