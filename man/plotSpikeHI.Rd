\name{plotSpikeHI}
\alias{plotSpikeHI}
\title{ plot hybridized target vs intensity}
\description{
  With final parameter setting, plot the amount of hybridized targets and intensities of
  calibration controls.  
}

\usage{
plotSpikeHI(spike, parameter, array = 1, xlab = "log(Hybridized)",
            ylab = "log(Intensity)", main = colnames(spike$R)[array],
            xlim = NULL, ylim = NULL, pch = 19, cex = 0.2, 
            cy5col = "black", cy3col = "black", noerror = TRUE,
            noepch = 19, noecex = 0.1, noecy5col = "lightpink",
            noecy3col = "lightblue", curve = TRUE, lwd = 1.5, 
            curvecy5col = "red", curvecy3col = "green", ...)
}

\arguments{
  \item{spike}{ a \code{\link[CALIB:spikelist]{SpikeList}} object.}

  \item{parameter}{ a \code{\link[CALIB:parameterlist]{ParameterList}} object.}

  \item{array}{ integer giving the array to be plotted. }

  \item{xlab}{ a title for the x axis.}

  \item{ylab}{ a title for the y axis.}

  \item{main}{ an overall title for the plot.}

  \item{xlim}{ the x limits (min,max) of the plot.}

  \item{ylim}{ the y limits of the plot.}

  \item{pch}{ an integer code for one of plotting characters or symbols
        for the spike data set. Default is 19.}

  \item{cex}{ a numerical value giving the amount by which the points which
        indicate spike data set should be scaled relative to the default. 
        Default is 0.2.}

  \item{cy5col}{ color of points for cy5. Default is black. }

  \item{cy3col}{ color of points for cy3. Default is black. }

  \item{noerror}{ a logical value. If it is \code{TRUE}, plot the amount of hybridized targets
        assuming equal spot capacities. Default is \code{TRUE}. }

  \item{noepch}{ pch for the points with equal spot capacities. Default is 19. }

  \item{noecex}{ cex for the points with equal spot capacities. Default is 0.1.}

  \item{noecy5col}{ color for the points with equal spot capacities of cy5. 
        Default is lightpink.}

  \item{noecy3col}{ color for the points with equal spot capacities of cy3.
        Default is lightblue.}

  \item{curve}{ a logical value. If it is \code{TRUE}, plot final parameter setting.
        Default is \code{TRUE}.}

  \item{lwd}{ width of the parameter curves. Default is 1.5.}

  \item{curvecy5col}{ color of the parameter curves for cy5. Default is red. }

  \item{curvecy3col}{ color of the parameter curves for cy3. Default is green. }

  \item{\dots}{ other graphical parameters can be used in function '\code{\link[graphics]{plot}}'. }
}

\details{
   The function plots hybridized targets vs measured intensities of one array. The argument 
   \code{array} gives the array index to be plotted. The function accepts the spike concentrations 
   from the argument \code{spike} and the estimated spot error for each spot from the argument 
   \code{parameter}. The hyrbidized targets for each spot can be calculated by the following 
   formula:
   formula.

   The argument \code{noerror} says whether or not the hybridized targets, which are calculated
   by the above mentioned formula assuming equal spot capacities, are plotted. If they are
   plotted, other arguments like \code{noepch},\code{noecex},\code{noecy5col} and \code{noecy3col} are used to 
   specify the type, the size and the color of the points.

   Estimated parameter curves can be shown on the plot. Since model parameters are different
   for two colors, two parameter curves are expected for one array. The funciton accepts 
   parameters of both colors from the argument \code{parameter}. If the curves are plotted, the 
   arguments \code{lwd}, \code{curvecy3col} and \code{curvecy5col} are used to specify the width and color
   of the curves.

   Details for the graphical parameters can be seen in function \code{\link[graphics]{plot}},
   \code{\link[graphics]{points}} and \code{\link[graphics]{curve}}.    
  
}

\value{
   A plot is created on the current graphics device.
}

\references{ Engelen, K., Naudts, B., DeMoor, B., Marchal, K. (2006) 
   A calibration method for estimating absolute expression levels from 
   microarray data. Bioinformatics  22: 1251-1258.}

\author{ Hui Zhao}

\examples{
# load data: spike and parameter
data(spike)
data(parameter)

# specify the array to be plotted.
array <- 1

# use the default values for other parameters
plotSpikeHI(spike,parameter,array=array)

}

\keyword{ hplot }
