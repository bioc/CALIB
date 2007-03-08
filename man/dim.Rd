\name{dim}
\alias{dim.RGList_CALIB}
\alias{dim.SpikeList}
\alias{length.RGList_CALIB}
\alias{length.SpikeList}
\title{Retrieve the Dimensions of an RGList\_CALIB or SpikeList object }

\description{
 Retrieve the number of rows (genes) and columns (arrays) for an 
 RGList\_CALIB or SpikeList
}

\usage{
 \method{dim}{RGList_CALIB}(x)
 \method{length}{RGList_CALIB}(x)
}

\arguments{

  \item{x}{ an object of class \code{RGList_CALIB} or \code{SpikeList}.}
}

\details{
  Microarray data objects share many analogies with ordinary
  matrices in which the rows correspond to spots or genes and the
  columns to arrays. These methods allow one to extract the size of
  microarray data objects in the same way that one would do for
  ordinary matrices.

  A consequence is that row and column commands \code{\link[base:nrow]{nrow(x)}},
  \code{\link[base:ncol]{ncol(x)}} and so on also work.
}

\value{
  Numeric vector of length 2. The first element is the number of rows(genes)
  and the second is the number of columns(arrays). 
}

\references{ \code{\link[limma]{dim}} in limma package}

\author{ Hui Zhao}

\seealso{ 
  \code{\link[base]{dim}} in the base package

  \code{\link[limma]{dim}} in the limma package
}

\examples{

# for RGList_CALIB
R <- G <- matrix(1:8,4,2)
rownames(R) <- rownames(G) <- c("g1","g2","g3","g4")
colnames(R) <- colnames(G) <- c("a1","a2")
RG <- new("RGList_CALIB",list(R=R,G=G))

dim(RG)

# for SpikeList
SR <- SG <- matrix(1:8,4,2)
rownames(SR) <- rownames(SG) <- c("s1","s2","s3","s4")
colnames(SR) <- colnames(SG) <- c("a1","a2")
spike <- new("SpikeList",list(R=SR,G=SG))

dim(spike)
}

\keyword{ array }
