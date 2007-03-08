\name{estimateParameter}
\alias{estimateParameter}
\title{Estimate model parameter from spikes }

\description{
 Estimate the calibration model parameters according to the known concentration and 
 the measured intensities of external control spikes on each array. 
}

\usage{
 estimateParameter(spike, RG, bc = FALSE, area = TRUE, errormodel = "M")
}

\arguments{
  \item{spike}{ a \code{SpikeList} object.}

  \item{RG}{ a \code{RGList_CALIB} object.}

  \item{bc}{ a logical value. \code{TRUE} means background corrected measured 
        intensities are used. Default is \code{FALSE}.}

  \item{area}{ a logical value. \code{TRUE} means spot area is used to calculate 
        measured intensities. Namly, measured intensities are calculated by 
        foreground intensities(or background corrected intensities, if bc is 
        \code{TRUE} ) multiply spot area. \code{FALSE} means spot area is not used.
        Default is \code{TRUE} .}

  \item{errormodel}{ a character to indicate the distribution of spot capacity. 
        "A" means spot capacity is additive. "M" means spot capacity is multiplicative.
        Default is "M". }
}

\details{
  This function estimates calibration model parameters. In this function, 
  the model parameters are estimated separately for each microarray, 
  based on the measured intensities of the external control spikes and 
  their known concentration in the hybridization solution. It accepts 
  spike measured intensities and concentration from \code{spike} argument, 
  which is an object of \code{\link[CALIB:spikelist]{SpikeList}} class. 

  It supports different ways to calculate the measured intensities. 
  Arguments \code{bc} and \code{area} are logical and their combinations is used 
  for specifying four differents ways. \code{bc} indicates using background correction
  or not. \code{area} indicates multiplying spot area or not. The default value
  of these two arguments are \code{bc} = FALSE and \code{area} = TRUE. 

  The argument \code{errormodel} is to specify the distribution of spot capacity 
  of each array. The spot capacity is either additive or multiplicative.
  Whichever distribution is more appropriate will depend largely on the type of 
  microarray slide and spotting procedure used. The spot parameters mus and sigmas
  can be considered equal for all measurements of a single array.

  The argument \code{RG} is for calculating the maximum intensity of each array. These 
  maximum intensities are used to estimate the upper saturation level of each array.

  More details please refer to the reference literature.
}

\value{
  An \code{ParameterList} object containing the components:
  
  \item{MuS}{ matrix containing MuS for each array.}

  \item{Ka}{ matrix containing Ka for each array.}

  \item{P1}{ matrix containing P1 of each dye for each array.  }

  \item{P2}{ matrix containing P2 of each dye for each array.}

  \item{SigmaA}{ matrix containing sigma additive for each array.}

  \item{SigmaM}{ matrix containing sigma multiplicative for each array.}

  \item{SigmaS}{ matrix containing sigma spoterror for each array.}

  \item{SpotError}{ matrix containing the spot error of each spot for each array. }

  \item{Method}{ boolean values indicating the way to calculate the measured intensities.}

  \item{ErrorModel}{ character \code{"M"} or \code{"A"} to indicate the type of spot 
        capacity distribution.}
}

\references{ Engelen, K., Naudts, B., DeMoor, B., Marchal, K. (2006) 
   A calibration method for estimating absolute expression levels from 
   microarray data. Bioinformatics  22: 1251-1258.}

\author{ Hui Zhao}

\examples{

# load data: RG and spike
data(RG)
data(spike)

# for the measured itensities, take the default bc=FALSE and area=TRUE.
# use multiplicative spot error model
parameter <- estimateParameter(spike,RG)

}

\keyword{ optimize }% at least one, from doc/KEYWORDS
