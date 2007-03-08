\name{plotSpikeRG}
\alias{plotSpikeRG}
\title{ plot spike intensity R vs G }

\description{
  plot red intensity vs green intensity of spikes.
}

\usage{
plotSpikeRG(spike,parameter,RG,array = 1, bc = FALSE, area = TRUE,
            xlab = "log(Rintensity)", ylab = "log(Gintensity)",
            main = colnames(spike$R)[array], onlycalib = FALSE, 
            xlim = NULL, ylim = NULL, pch = 19, cex = 0.3, col = "black", 
            allpch = 19, allcex = 0.05, allcol = "lightgrey", diag = TRUE,
            diagcol = "grey", diaglwd = 1, curvecol = "blue",
            curvelwd = 1.5, calibtype = 1, adjusttype = 4, ...)
}

\arguments{
  \item{spike}{ a \code{\link[CALIB:spikelist]{SpikeList}} object.}
  
  \item{parameter}{a \code{\link[CALIB:parameterlist]{ParameterList} object.}.
       If \code{parameter} argument is sepcified, model curves are shown on the plot. }

  \item{RG}{a \code{\link[CALIB:RGList_CALIB]{RGList_CALIB} object}. If \code{parameter}
       argument is specified, this argument is obligated. More description in Detail section.}

  \item{array}{ integer giving the array to be plotted. }

  \item{bc}{  a logical value. \code{TRUE} means background corrected measured 
        intensities are used. Default is \code{FALSE}.}

  \item{area}{  a logical value. \code{TRUE} means spot area is used to calculate 
        measured intensities. Namly, measured intensities are calculated by 
        foreground intensities(or background corrected intensities, if bc is 
        \code{TRUE}) multiply spot area. \code{FALSE} means spot area is not used.
        Default is \code{TRUE}.}

  \item{xlab}{ a title for the x axis.}

  \item{ylab}{ a title for the y axis.}

  \item{main}{ an overall title for the plot.}

  \item{onlycalib}{ a logical value. \code{TRUE} means only the calibration controls are on 
        the plot. \code{FALSE} means to plot all the spikes}

  \item{xlim}{ the x limits (min,max) of the plot.}

  \item{ylim}{ the y limits (min,max) of the plot.}

  \item{pch}{ an integer code for one of a set of plotting characters
        or symbols for the spike data set. Default is 21.}

  \item{cex}{ a numerical value giving the amount by which the points which 
        indicate spike data set should be scaled relative to the default.
        Default is 0.3.}

  \item{col}{ the color of the points indicating spike data set. Default is black.}
 
  \item{allpch}{ an integer code for one of a set of plotting characters
        or symbols for the spike data set. Default is 19. }

  \item{allcex}{ a numerical value giving the amount by which the points which 
        indicate all data set should be scaled relative to the default. Default is 0.05. }

  \item{allcol}{the color of the points indicating all data set. Default is lightgrey.}

  \item{diag}{ a logical value. Add diagonal on the plot if it is \code{TRUE}.
        Default is \code{TRUE}.}

  \item{diagcol}{ the color of the diagonal. Default is grey.}

  \item{diaglwd}{ the width of the diagonal. Default is 1.}
 
  \item{curvecol}{ the color of the model curves. Default is blue.}

  \item{curvelwd}{the width of the model curves for calibration
        control spikes. Default is 1.5.}
 
  \item{calibtype}{the line type of the model curves for 
        calibration control spikes. Default is 1.}

  \item{adjusttype}{the line type of the model curves (using parameter after
       adjustment) for calibration control spikes. Default is 4.}

  \item{\dots}{ other graphical parameters can be used in function \code{\link[graphics]{plot}}. }
}

\details{
  The function plots red vs green measured intensities of spikes of one array. 
  The argument \code{array} gives the array index to be plotted. 

  If parameter estimation is done, the model curves can be plotted by giving
  the argument \code{parameter}. And if the argument \code{parameter} is specified,
  the argument \code{RG} is obligated. The combination of these two arguments is
  used to compare how the model fits to the spike and to the whole data set. More
  details about the usage of this function refers to the readme file of this package.

  It supports different ways to calculate the measured intensities. 
  Arguments \code{bc} and \code{area} are logical and their combinations are used 
  for specifying four differents ways. \code{bc} indicates using background correction
  or not. \code{area} indicates multipling spot area or not. The default value
  of these two arguments are \code{bc} = FALSE and \code{area} = TRUE. 

  see other graphic functions for the other arguments. 
}

\value{
  A plot is created on the current graphics device.
}

\author{ Hui Zhao}

\examples{
# load data: spike
data(spike)

# specify the array to be plotted.
array <- 1

# use the default values for other arguments
plotSpikeRG(spike,array=array)

# after parameter estimation, the model curves can be shown on the plot.
data(parameter)
data(RG)
plotSpikeRG(spike,parameter,RG)

}

\keyword{ hplot }