\name{merge}
\alias{merge.RGList_CALIB}
\alias{merge.SpikeList}
\title{Merge RGList\_CALIB or SpikeList objects}

\description{
 Merge two \code{RGList_CALIB} objects or two \code{SpikeList} objects 
 in possibly irregular order. 
}

\usage{
 \method{merge}{RGList_CALIB}(x,y,...)
}

\arguments{

  \item{x}{ an \code{RGList_CALIB} object or an \code{SpikeList} object.}

  \item{y}{ corresponding \code{RGList_CALIB} object or \code{SpikeList} 
        object. Has the same genes or spikes as \code{x}, possibly in a different
        order,but with different arrays. }

  \item{\dots}{ other arguments can be used in \code{\link[base]{merge}} in the base packge.}

}

\details{
  \code{\link[CALIB:rglist_CALIB]{RGList_CALIB}} and \code{\link[CALIB:spikelist]{SpikeList}} 
  objects are list objects containing numeric matrices with the same dimensions. 
  The RGLists\_CALIB or SpikeLists are merged by merging each of the components by row 
  names or, if there are no row names, by IDs in the \code{genes} component. Unlike
  when using \code{cbind}, row names are not required to be in the same order
  or to be unique. In the case of repeated row names, the order of the rows with 
  repeated names in preserved. This means that the first occurrence of each name 
  in \code{x$R} is matched with the first occurrence of the same name in \code{y$R},
  the second with the second, and so on. The final vector of row names is the same 
  as in \code{x}.
}

\note{
  If the \code{RGList_CALIB} or \code{SpikeList} objects contain the same number of genes 
  or spikes in the same order then the appropriate function to combine them is 
  \code{\link[CALIB:cbind]{cbind}} rather than \code{merge}. 
}

\value{
  An merged object of the same class as \code{x} and \code{y} with the same components as
  \code{x}. Components matrices have the same row names as in \code{x} but columns from
  \code{y} as well as \code{x}.
}

\references{ \code{\link[limma]{merge}} in limma package}

\author{ Hui Zhao}

\seealso{ 
  \code{\link[base]{merge}} in the base package

  \code{\link[limma]{merge}} in the limma package
}

\examples{

R1 <- G1 <- matrix(1:8,4,2)
rownames(R1) <- rownames(G1) <- c("g1","g1","g2","g3")
colnames(R1) <- colnames(G1) <- c("a1","a2")
RG1 <- new("RGList_CALIB",list(R=R1,G=G1))

R2 <- G2 <- matrix(9:16,4,2)
rownames(R2) <- rownames(G2) <- c("g2","g3","g1","g1")
colnames(R2) <- colnames(G2) <- c("a3","a4")
RG2 <- new("RGList_CALIB",list(R=R2,G=G2))

RG12 <- merge(RG1,RG2)
RG21 <- merge(RG2,RG1)

}

\keyword{ manip }
