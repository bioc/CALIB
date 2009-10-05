\name{subsetting}
\alias{[.RGList_CALIB}
\alias{[.SpikeList}
\alias{[.ParameterList}
\title{Subset of RGList\_CALIB, SpikeList or ParameterList object}

\description{
 Exact a subset of an RGList\_CALIB, SpikeList or ParameterList object.
}

\usage{
 \method{[}{RGList_CALIB}(object, i, j, \ldots)
}

\arguments{

  \item{object}{ an object of class \code{RGList_CALIB} or \code{SpikeList}.}
  \item{i,j}{subscripts. \code{i} is the subscripts of the spots 
             and \code{j} is the subscripts of the arrays.}
  \item{\ldots}{not used}

}

\details{
  \code{i},\code{j} may take any values acceptable for the matrix
components of \code{object}. 
}

\value{
 An object of the same class as \code{object} containing the data with specified
 subset of spots and arrays.
}

\references{ \code{\link[limma]{subsetting}} in limma package}

\author{ Hui Zhao}

\seealso{ 
   \code{\link[limma]{subsetting}} in the limma package
}

\examples{

# for RGList_CALIB
R <- G <- matrix(1:8,4,2)
rownames(R) <- rownames(G) <- c("g1","g2","g3","g4")
colnames(R) <- colnames(G) <- c("a1","a2")
RG <- new("RGList_CALIB",list(R=R,G=G))

RG[1:2,]
RG[,1:2]
RG[1:2,1:2]
}

\keyword{ array }
