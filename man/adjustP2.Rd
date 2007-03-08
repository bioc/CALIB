\name{adjustP2}
\alias{adjustP2}
\title{ Adjust model parameter P2 }
\description{
Adjust the calibration model parameter P2 according to the measured intensities
of all clones spotted on the array.
}
\usage{
adjustP2(RG, parameter, arrayindex = arrayindex, colorindex = colorindex)
}
\arguments{
  \item{RG}{ an \code{RGList} object.}
  \item{parameter}{ a \code{ParameterList} object. }
  \item{arrayindex}{ integer vector specifying the index of the arrays of whose 
        the parameter P2 needed to be adjusted. }
  \item{colorindex}{ integer vector specifying which color needed to be adjusted.}
}
\details{
\code{RG} is an \code{\link[CALIB:rglist_CALIb]{RGList_CALIB}} object which contains all the
experimental data. 
\code{parameter} is the return result of function \code{\link[CALIB:estimateParameter]{estimateParameter}.}
\code{arrayindex} is an integer vector. It gives the index of the arrays whose P2 is 
needed to be adjusted.
\code{colorindex} is an integer vector. It gives the color needed to be adjusted.
\code{1} means red P2 and \code{2} means green P2. 
}
\value{
 It returns a \code{\link[CALIB:parameterlist]{ParameterList}} object with a adjusted
P2 compared to the input argument \code{parameter}. 
}
\note{
The user should decide on which array and which color the adjustment is needed.
Therefore it is important to specify the right array index and color index.
There is no check on this in the function. 
}
\author{ Hui Zhao}
\examples{
# load data: RG and parameter:
data(RG)
data(parameter)

# adjust P2
parameter_new <- adjustP2(RG,parameter,arrayindex=c(1,2),colorindex=c(2,2))
}
\keyword{ optimize }