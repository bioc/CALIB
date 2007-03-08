\name{plotSpikeSpotError}
\alias{plotSpikeSpotError}
\title{ plot spot error of spikes.}

\description{
  Plot spot error of spikes. Three types of plots are avaible: histogram, 
  boxplot and density. 
}

\usage{
plotSpikeSpotError(parameter, array = 1, plottype = "hist", width = 1, 
                   plotnames = NULL, main = NULL, ...)
}

\arguments{
  \item{parameter}{ a \code{\link[CALIB:parameterlist]{ParameterList}} object.}

  \item{array}{integer giving the array to be plotted. }

  \item{plottype}{ string giving the type of plot. }

  \item{width}{ needed for density plot. This exists for compatibility with S; 
        if given, and \code{bw} is not, will set \code{bw} to \code{width}
        if this is a character string, or to a kernel-dependent multiple of 
        \code{width} if this is numeric. Default is 1.}

  \item{plotnames}{ needed for boxplot. group labels which will be printed 
        under each boxplot. }

  \item{main}{ an overall title for the plot.}

  \item{\dots}{ other parameters can be used according to the plottype user specified.}
}

\details{
  The function plots spot error of one array on different types of plots. Three types, 
  which are histogram, boxplot and density function, are available now. The argument 
  \code{plottype} is used for giving the plot type. It should be one of the following three 
  types: "hist", "boxplot" and "dens". The argument \code{array} gives the array index to 
  be plotted. 

  The function accepts estimated spot error from the argument \code{parameter}.   
}

\value{
  A plot is created on the current graphics device.
}

\author{ Hui Zhao}

\seealso{ 
\code{\link[graphics]{hist}}, \code{\link[graphics]{boxplot}}
and \code{\link[graphics]{plot}} in the graphics package.

\code{\link[stat]{density}} in the stat package.}

\examples{
# load data: parameter
data(parameter)

# specify the array to be plotted.
array <- 1

# plot histogram
plotSpikeSpotError(parameter,array=array,plottype="hist")
# plot boxplot
plotSpikeSpotError(parameter,array=array,plottype="boxplot",plotnames=NULL)
# plot density function
plotSpikeSpotError(parameter,array=array,plottype="dens",width=1)

}

\keyword{ hplot }
