\name{dimnames}
\alias{dimnames.RGList_CALIB}
\alias{dimnames.SpikeList}
\title{Retrieve the Dimension Names of an RGList\_CALIB or SpikeList object.} 

\description{
 Retrieve the the dimension names of an RGList\_CALIB object or an SpikeList object
}

\usage{
 \method{dimnames}{RGList_CALIB}(x)
}

\arguments{

  \item{x}{ an object of class \code{RGList_CALIB} or \code{SpikeList}.}
}

\details{
  The dimension names of a micrroarry object or a spike object are the same as
  those of the most important matrix component of that object.

  A consequence is that row and column command \code{rownames} and \code{colnames} also work.
}

\value{
  Either \code{NULL} or a list of length 2. If the value is a list, its componets are either
  \code{NULL} or a character vector with the length of the appropriate dimension of \code{x}. 
  If the list component is not \code{NULL}, the first field of the list indicates rownames 
  and the second field indicates colnames.
}

\references{ \code{\link[limma]{dimnames}} in limma package}

\author{ Hui Zhao}

\seealso{ 
  \code{\link[base]{dimnames}} in the base package

  \code{\link[limma]{dimnames}} in the limma package
}

\examples{

# for RGList_CALIB
R <- G <- matrix(1:8,4,2)
rownames(R) <- rownames(G) <- c("g1","g2","g3","g4")
colnames(R) <- colnames(G) <- c("a1","a2")
RG <- new("RGList_CALIB",list(R=R,G=G))

dimnames(RG)

# for SpikeList
SR <- SG <- matrix(1:8,4,2)
rownames(SR) <- rownames(SG) <- c("s1","s2","s3","s4")
colnames(SR) <- colnames(SG) <- c("a1","a2")
spike <- new("SpikeList",list(R=SR,G=SG))

dimnames(spike)
}

\keyword{ array }
