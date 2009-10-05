\name{parameter}
\alias{parameter}
\docType{data}

\title{ Calibration Model parameter: ParameterList Example }

\description{
This \code{\link[CALIB:ParameterList-class]{ParameterList}} object
represents calibration  
model parameters of two arrays. It was estimated from the dataset
\code{\link{spike}} 
by the function \code{\link{estimateParameter}} in the CALIB package. 
}

\usage{data(parameter)}

\format{
\code{parameter} is a \code{ParameterList} object containing the following list 
components:\$MuS,\$Ka,\cr\$P1,\$P2,\$SigmaA,\$SigmaM,\$SigmaS,\$SpotError,\$Method
with two subfields \$Method\$BC and \$Method\$Area, \$ErrorModel and \$genes.
Among these parameters,P1,P2 and SigamA are different for different dyes.
}

\references{
dataset \code{spike}.

Engelen, K., Naudts, B., DeMoor, B., Marchal, K. (2006) 
A calibration method for estimating absolute expression levels from 
microarray data. Bioinformatics  22: 1251-1258.

Hilson,P.,et al. (2004) Versatile gene-specific sequence tags for 
Arabidopsis functional genomics: transcript profiling and reverse
genetics applications. Genome Res. 14, 2176-2189.

}

\examples{
data(parameter)
}

\keyword{datasets}
