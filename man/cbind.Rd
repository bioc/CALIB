\name{cbind}
\alias{cbind.RGList_CALIB}
\alias{cbind.SpikeList}
\alias{cbind.ParameterList}
\alias{rbind.RGList_CALIB}
\alias{rbind.SpikeList}
\title{Combine RGList\_CALIB, SpikeList or ParameterList objects}
\description{
 Combine a series of \code{RGList\_CALIB} objects or a series of \code{SpikeList} 
 objects or a series of\cr \code{ParameterList} objects. 
}
\usage{
 \method{cbind}{RGList_CALIB}(..., deparse.level = 1)
 \method{rbind}{RGList_CALIB}(..., deparse.level = 1) 
}
\arguments{

  \item{\dots}{ \code{RGList_CALIB} objects, \code{SpikeList} objects
        or \code{ParameterList} objects .}

  \item{deparse.level}{see \code{\link[base]{cbind}} in base package.}
}
\details{
  \code{cbind} combines data objects assuming the same gene lists but
  different arrays. \code{rbind} combines data objects assuming
  equivalent arrays, i.e., the same RNA targets, but different
  genes.

  For \code{ParameterList} objects, only \code{cbind} is available, because it makes 
  no sense to rbind parameter.

  For \code{RGList_CALIB} objects and \code{SpikeList} objects, cbind and rbind are both 
  available. For \code{cbind}, the matrices of expression data from the 
  individual objects are cbinded. The data.frames of target information, 
  if they exist, are rbinded. The combined data object will preserve
  any additional components or attributes found in the first object
  to be combined. For \code{rbind}, the matrices of expression data are
  rbinded while the target information, in any, is unchanged.
}
\value{
  An \code{\link[CALIB:RGList_CALIB-class]{RGList_CAILB}}, a
  \code{\link[CALIB:SpikeList-class]{SpikeList}} or  
  a \code{\link[CALIB:ParameterList-class]{ParameterList}} object
  holding data from all the arrays and all genes from the individual objects.
}
\references{ \code{\link[limma]{cbind}} in limma package}
\author{ Hui Zhao}
\seealso{ 
  \code{\link[base]{cbind}} in the base package

  \code{\link[limma]{cbind}} in the limma package
}
\examples{

R1 <- G1 <- matrix(1:8,4,2)
rownames(R1) <- rownames(G1) <- c("g1","g2","g3","g4")
colnames(R1) <- colnames(G1) <- c("a1","a2")
RG1 <- new("RGList_CALIB",list(R=R1,G=G1))

R2 <- G2 <- matrix(9:16,4,2)
rownames(R2) <- rownames(G2) <- c("g1","g2","g3","g4")
colnames(R2) <- colnames(G2) <- c("a3","a4")
RG2 <- new("RGList_CALIB",list(R=R2,G=G2))

RG <- cbind(RG1,RG2)

}
\keyword{ manip }
